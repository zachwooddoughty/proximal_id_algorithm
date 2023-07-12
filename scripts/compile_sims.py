import argparse
import csv
import glob
import json
import logging
import os
import re

from collections import defaultdict

import numpy as np

defaults = dict(
    z_bin=0, m_bin=0,
    ay_effect="NA", zmw_effect="NA",
    uw_effect="NA", uz_effect="NA",
    n_bootstrap=1,
)

methods = ["backdoor", "naive", "simple", "pfd"]
raw_metrics = ["estimate", "width", "coverage"]

method_names = {"backdoor": "Oracle Backdoor", "naive": "Naive Front-Door",
                "simple": "Simple Proximal", "pfd": "Proximal Front-Door"}
column_names = {"n_samples": "Sample Size", "ay_effect": r"$A \to Y$",
                "zmw_effect": r"$Z \to M \to W$", "u_effect": r"$U \to \{W,Z\}$"}

# Example filename
# z_bi0-m_bi0-ay_e0-uz_eNA-uw_eNA-n_bo64-dgp1-data1-samples16000.csv


def get_default_dict(depth=1):
  assert type(depth) == int
  if depth <= 1:
    return defaultdict(list)
  else:
    return defaultdict(lambda: get_default_dict(depth - 1))


def metric_agg_name(metric, agg):
  if metric == "pab.estimate":
    if agg in ["mean", "meanabs"]:
      return "Percent Absolute Bias"
  if metric == "estimate":
    if agg == "mean":
      return "Mean Bias"
    if agg == "meanabs":
      return "Mean Absolute Bias"
  if metric == "coverage":
    return "Bootstrap Interval Coverage"
  if metric == "width":
    return "Bootstrap Interval Width"

  return "{}.{}".format(metric, agg)


def print_table(data, column, metric, headers, **kwargs):

  true_effects = kwargs.get("true_effects", None)
  agg = kwargs.get("agg", "mean")
  latex = kwargs.get("latex", False)
  prec = kwargs.get("prec", 3)
  cell_width = kwargs.get("cell_width", 7)

  if latex:
    method_width = max(20, max(map(len, method_names.values())))
    row = ["{:{width}s}".format(column_names[column], width=method_width)]
  else:
    title = "{} {} {}".format(column, metric, agg)
    method_width = max(20, max(map(len, method_names.values())), len(title))
    row = ["{:{width}s}".format(title, width=method_width)]

  for header in headers:
    row.append("{:{width}}".format(header, width=cell_width))
  if latex:
    print(r"{:s} & \multicolumn{{{:d}}}{{c}}{{{:s}}} \\".format(
        "Metric", len(headers), metric_agg_name(metric, agg)))
    print(" & ".join(row) + r" \\")
    print(r"\midrule")
  else:
    print(" ".join(row))

  if agg == "mean":
    def agg_func(data):
      dgp_means = []
      for dgp in data:
        mean = np.mean(data[dgp])
        if "pab" in metric:
          assert true_effects is not None, "pab must pass true_effect"
          mean /= true_effects[dgp]
        dgp_means.append(mean)
      return np.mean(dgp_means)

  elif agg == 'std':
    def agg_func(data):
      return np.mean([np.std(vals) for dgp, vals in data.items()])
  elif agg in "meanabs":
    def agg_func(data):
      dgp_means = []
      for dgp in data:
        mean = np.abs(np.mean(data[dgp]))
        if "pab" in metric:
          assert true_effects is not None, "pab must pass true_effect"
          mean /= np.abs(true_effects[dgp])
        dgp_means.append(mean)
      return np.mean(dgp_means)
  else:
    raise ValueError("Unknown aggregation function {}".format(agg))

  for method in methods:
    row = ["{:{width}s}".format(method_names[method], width=method_width)]
    for header in headers:
      cell_data = data[method][header]
      val = agg_func(cell_data)
      row.append("{:{width}.{prec}f}".format(val, width=cell_width, prec=prec))
    if latex:
      print(" & ".join(row) + r" \\")
    else:
      print(" ".join(row))

  if latex:
    print(r"\bottomrule")


def build_table(indir, columns, aggs="mean,std", filters={}, latex=False, prec=3):

  table = get_default_dict(4)
  counts = defaultdict(list)
  all_dgps = set()
  all_datasets = set()

  assert len(columns) == 1
  errors = 0
  not_done = 0

  metrics = raw_metrics
  if filters.get("n_bootstrap", 1) <= 1:
    metrics = ["estimate"]

  true_effects = {}
  dgp_regex = re.compile(r"dgp(\d+)")
  data_regex = re.compile(r"data(\d+)")
  column, column_values = next(iter(columns.items()))
  for column_value in column_values:
    fns = match_job_params(indir, **{column: column_value}, **filters)
    for fn in fns:
      try:
        results = read_job(fn)
      except ValueError as e:
        errors += 1
        logging.warning("{} on {}".format(e, fn))
        continue
      except IndexError as e:
        errors += 1
        logging.warning("{} on {}".format(e, fn))
        continue
      except StopIteration:
        not_done += 1
        logging.warning("StopIter on {}".format(fn))
        continue

      dgp_i = int(dgp_regex.search(fn).group(1))
      data_i = int(data_regex.search(fn).group(1))
      all_dgps.add(dgp_i)
      all_datasets.add(data_i)
      true_effect = results["dgp_params"]["true_effect"][0]
      true_effects[int(dgp_i)] = true_effect

      for method in methods:
        for metric in metrics:
          val = results[method][metric]
          table[metric][method][column_value][dgp_i].append(val)
          counts["{}_{}_{}".format(metric, method, column_value)].append(
              (dgp_i, data_i))

          if metric in ["estimate", "mean"]:
            pab_metric = "pab.{}".format(metric)
            table[pab_metric][method][column_value][dgp_i].append(val)

  count_list = [len(x) for x in counts.values()]
  if len(count_list) == 0:
    print("No results")
    return
  print("{} errors; {} not done; mean {:.1f} max {:.3f}".format(
      errors, not_done,
      np.mean(count_list), np.max(count_list)))

  all_dgps = list(range(1, 5))
  all_datasets = list(range(1, 65))

  expected = len(all_dgps) * len(all_datasets)
  missing = defaultdict(int)
  if np.mean(count_list) != expected:
    for key in counts:
      if len(counts[key]) != expected:
        for dgp_i in all_dgps:
          for data_i in all_datasets:
            tuple_key = (dgp_i, data_i)
            if tuple_key not in counts[key]:
              metric, method, column_value = key.split("_")
              missing['"{}-{}-{}"'.format(column_value, *tuple_key)] += 1

    if len(missing) > 0:
      print(" ".join(list(missing.keys())))

  for metric in sorted(table):
    for agg in aggs.split(","):
      kwargs = dict(agg=agg, latex=latex, prec=prec)
      if 'pab' in metric:
        kwargs["true_effects"] = true_effects
      print_table(table[metric], column, metric,
                  sorted(column_values), **kwargs)


def read_job(fn):
  d = {}
  with open(fn) as inf:
    d["job_params"] = json.loads(inf.readline()[3:])
    d["dgp_params"] = json.loads(inf.readline()[3:])
    reader = csv.reader(inf)
    headers = next(reader)
    data = [row for row in reader]

  for row in data:
    method = row[0]
    d[method] = {}
    for i in range(1, len(headers)):
      header = headers[i]
      d[method][header] = float(row[i])

  return d


def match_job_params(indir, **kwargs):
  keys = {}
  for key in defaults:
    keys[key] = kwargs.get(key, defaults[key])

  if "u_effect" in kwargs:
    keys["uw_effect"] = kwargs["u_effect"]
    keys["uz_effect"] = kwargs["u_effect"]

  fn_params = []
  for key in ["z_bin", "m_bin", "ay_effect",
              "uz_effect", "uw_effect", "zmw_effect",
              "n_bootstrap"]:
    val = keys[key]

    fn_params.append("{}{}".format(key[:4], val))

  fn_template = "-".join(fn_params) + "-dgp*-data*-samples{samples}.csv"
  if "n_samples" in kwargs:
    fn_template = fn_template.format(samples=kwargs["n_samples"])
  else:
    fn_template = fn_template.format(samples="*")

  for fn in glob.glob(os.path.join(indir, fn_template)):
    yield fn

  # if zmw is None, also iterate through examples with zmw_NA
  if keys["zmw_effect"] == "NA":
    second_template = fn_template.replace("zmw_NA-", "")
    for fn in glob.glob(os.path.join(indir, second_template)):
      yield fn


def main():
  parser = argparse.ArgumentParser()
  parser.add_argument("indir", type=str)
  parser.add_argument("table", type=int, choices=[1, 2, 3, 4])
  parser.add_argument("--n_samples", type=int, default=4000)
  parser.add_argument("--n_bootstrap", type=int, default=64)
  parser.add_argument("--zm_bin", action='store_true')
  parser.add_argument("--dgp", type=int, default=None)
  parser.add_argument("--agg", type=str, default="mean,std")
  parser.add_argument("--latex", action='store_true')
  parser.add_argument("--prec", type=int, default=3)
  args = parser.parse_args()

  filters = {
      "n_samples": args.n_samples, "n_bootstrap": args.n_bootstrap,
      "z_bin": int(args.zm_bin), "m_bin": int(args.zm_bin),
  }
  if args.dgp is not None:
    filters["dgp"] = args.dgp

  if args.table == 1:
    filters.pop("n_samples")
    build_table(args.indir, {"n_samples": [4000, 16000, 64000]},
                filters=filters, aggs=args.agg, latex=args.latex, prec=args.prec)
  elif args.table == 2:
    build_table(args.indir, {"ay_effect": [0, 0.2, 0.4, 0.8]},
                filters=filters, aggs=args.agg, latex=args.latex, prec=args.prec)
  elif args.table == 3:
    build_table(args.indir, {"zmw_effect": [0, 0.2, 0.4, 0.8]},
                filters=filters, aggs=args.agg, latex=args.latex, prec=args.prec)
  elif args.table == 4:
    build_table(args.indir, {"u_effect": [0, 0.2, 0.4, 0.8]},
                filters=filters, aggs=args.agg, latex=args.latex, prec=args.prec)


if __name__ == "__main__":
  main()
