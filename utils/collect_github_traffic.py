#!/usr/bin/env python3
"""
collect_github_traffic.py — keep a permanent history of GitHub traffic.

GitHub's Insights -> Traffic only keeps 14 days and offers no download counter for a
repository. This script snapshots what GitHub does expose, merges it into CSV files that
grow over time, and adds the counters that *are* cumulative (Docker Hub pulls, release
asset downloads, stars/forks).

    python3 utils/collect_github_traffic.py                 # collect and merge
    python3 utils/collect_github_traffic.py --summary       # print what is stored, no network

Needs the `gh` CLI authenticated with push access to the repository (traffic endpoints
are not public). Run it at least every 14 days or days fall off the end and are lost.

Files, under metrics/ by default:
    traffic_clones.csv      date, clones, uniques
    traffic_views.csv       date, views, uniques
    referrers.csv           collected_on, referrer, views, uniques     (snapshot, 14-day window)
    paths.csv               collected_on, path, views, uniques         (snapshot, 14-day window)
    totals.csv              collected_on, docker_pulls, stars, forks, watchers, release_downloads

Daily rows are merged on the date: a later collection overwrites an earlier one for the
same day, because the day is only final once it has passed. Snapshot files append.
"""
import argparse
import csv
import datetime as dt
import json
import os
import subprocess
import sys
import urllib.request

REPO = "WeTGI-colab/SMART"
DOCKER = "monkiky/smart"
HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(HERE, "..", "metrics")


def gh(path):
    """One GitHub API call through the gh CLI, so no token is stored anywhere."""
    try:
        out = subprocess.run(["gh", "api", path], capture_output=True, text=True, check=True)
    except FileNotFoundError:
        sys.exit("gh not found: install the GitHub CLI (brew install gh)")
    except subprocess.CalledProcessError as e:
        msg = (e.stderr or "").strip().splitlines()
        sys.exit(f"gh api {path} failed: {msg[0] if msg else e}")
    return json.loads(out.stdout)


def docker_pulls(repo):
    try:
        with urllib.request.urlopen(f"https://hub.docker.com/v2/repositories/{repo}/",
                                    timeout=30) as fh:
            return json.load(fh).get("pull_count")
    except Exception as e:                       # the history is still worth keeping without it
        print(f"  note: Docker Hub unavailable ({e.__class__.__name__}), pulls left blank")
        return ""


def merge_daily(path, rows, columns):
    """Rows keyed by date; a later collection replaces an earlier one for the same day."""
    stored = {}
    if os.path.exists(path):
        with open(path, newline="") as fh:
            stored = {r["date"]: r for r in csv.DictReader(fh)}
    added = sum(1 for r in rows if r["date"] not in stored)
    stored.update({r["date"]: r for r in rows})
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, columns)
        w.writeheader()
        for date in sorted(stored):
            w.writerow({c: stored[date].get(c, "") for c in columns})
    return added, len(stored)


def append(path, rows, columns):
    new = not os.path.exists(path)
    with open(path, "a", newline="") as fh:
        w = csv.DictWriter(fh, columns)
        if new:
            w.writeheader()
        w.writerows(rows)


def collect(out_dir, repo, docker):
    os.makedirs(out_dir, exist_ok=True)
    today = dt.date.today().isoformat()
    print(f"collecting {repo} on {today}")

    clones = gh(f"repos/{repo}/traffic/clones")
    rows = [{"date": c["timestamp"][:10], "clones": c["count"], "uniques": c["uniques"]}
            for c in clones.get("clones", [])]
    added, total = merge_daily(os.path.join(out_dir, "traffic_clones.csv"), rows,
                               ["date", "clones", "uniques"])
    print(f"  clones: {clones.get('count')} in the last 14 days; "
          f"{added} new day(s), {total} day(s) on record")

    views = gh(f"repos/{repo}/traffic/views")
    rows = [{"date": v["timestamp"][:10], "views": v["count"], "uniques": v["uniques"]}
            for v in views.get("views", [])]
    added, total = merge_daily(os.path.join(out_dir, "traffic_views.csv"), rows,
                               ["date", "views", "uniques"])
    print(f"  views:  {views.get('count')} in the last 14 days; "
          f"{added} new day(s), {total} day(s) on record")

    for name, endpoint, key in [("referrers", "popular/referrers", "referrer"),
                                ("paths", "popular/paths", "path")]:
        items = gh(f"repos/{repo}/traffic/{endpoint}")
        append(os.path.join(out_dir, f"{name}.csv"),
               [{"collected_on": today, key: i.get(key) or i.get("title", ""),
                 "views": i["count"], "uniques": i["uniques"]} for i in items],
               ["collected_on", key, "views", "uniques"])
        print(f"  {name}: {len(items)} recorded")

    meta = gh(f"repos/{repo}")
    releases = gh(f"repos/{repo}/releases")
    downloads = sum(a.get("download_count", 0) for r in releases for a in r.get("assets", []))
    append(os.path.join(out_dir, "totals.csv"), [{
        "collected_on": today, "docker_pulls": docker_pulls(docker),
        "stars": meta.get("stargazers_count"), "forks": meta.get("forks_count"),
        "watchers": meta.get("subscribers_count"), "release_downloads": downloads,
    }], ["collected_on", "docker_pulls", "stars", "forks", "watchers", "release_downloads"])
    print(f"  totals: {meta.get('stargazers_count')} stars, {meta.get('forks_count')} forks, "
          f"{downloads} release-asset downloads")
    print(f"written to {os.path.abspath(out_dir)}")


def summary(out_dir):
    for name, value_col in [("traffic_clones.csv", "clones"), ("traffic_views.csv", "views")]:
        path = os.path.join(out_dir, name)
        if not os.path.exists(path):
            print(f"{name}: nothing collected yet")
            continue
        with open(path, newline="") as fh:
            rows = list(csv.DictReader(fh))
        total = sum(int(r[value_col]) for r in rows if r[value_col])
        uniq = sum(int(r["uniques"]) for r in rows if r["uniques"])
        gaps = 0
        dates = sorted(dt.date.fromisoformat(r["date"]) for r in rows)
        for a, b in zip(dates, dates[1:]):
            gaps += (b - a).days - 1
        print(f"{value_col}: {total} over {len(rows)} day(s), {uniq} unique, "
              f"{dates[0]} to {dates[-1]}" + (f", {gaps} day(s) missing" if gaps else ""))
    path = os.path.join(out_dir, "totals.csv")
    if os.path.exists(path):
        with open(path, newline="") as fh:
            rows = list(csv.DictReader(fh))
        first, last = rows[0], rows[-1]
        print(f"docker pulls: {last['docker_pulls']} on {last['collected_on']}"
              + (f" (was {first['docker_pulls']} on {first['collected_on']})"
                 if len(rows) > 1 else ""))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=OUT_DIR, help="directory for the CSV files")
    ap.add_argument("--repo", default=REPO)
    ap.add_argument("--docker", default=DOCKER)
    ap.add_argument("--summary", action="store_true", help="print what is stored, no network")
    a = ap.parse_args()
    summary(a.out) if a.summary else collect(a.out, a.repo, a.docker)
