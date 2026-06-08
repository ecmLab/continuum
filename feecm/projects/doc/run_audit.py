#!/usr/bin/env python3
"""Regenerate project_audit.tex by running --check-input on every .i file
under projects/ecmTest and projects/tecmTest."""
import os, re, subprocess, sys
from concurrent.futures import ProcessPoolExecutor
from datetime import date

FEECM = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
APPS = {
    "ecmTest":  os.path.join(FEECM, "ecm_test", "ecm-opt"),
    "tecmTest": os.path.join(FEECM, "tecm_test", "tecm_test-opt"),
}

def discover():
    files = []
    for app in APPS:
        root = os.path.join(FEECM, "projects", app)
        for dirpath, dirnames, filenames in os.walk(root):
            dirnames[:] = [d for d in dirnames if d not in (".jitcache", ".claude", "build")]
            for fn in filenames:
                if fn.endswith(".i"):
                    files.append((app, os.path.join(dirpath, fn)))
    return sorted(files, key=lambda x: x[1])

def categorize(out):
    o = re.sub(r"\x1b\[[0-9;]*m", "", out)
    m = re.search(r"'([A-Za-z0-9_]+)' is not a registered object", o)
    if m: return ("Missing class", m.group(1))
    m = re.search(r"Unable to open file \"([^\"]+)\"", o)
    if m: return ("Missing file", os.path.basename(m.group(1)))
    m = re.search(r"unable to open (?:file )?'?([^\s']+\.(?:e|msh|csv|cpa\.gz|exo))'?", o, re.I)
    if m: return ("Missing file", os.path.basename(m.group(1)))
    if re.search(r"(does not exist|could not.*open|No such file).*\.(e|msh|csv|gz|exo)", o, re.I):
        m = re.search(r"([^\s\"']+\.(?:e|msh|csv|cpa\.gz|exo))", o)
        return ("Missing file", os.path.basename(m.group(1)) if m else "data file")
    if "checkpoint" in o.lower() and ("restart" in o.lower() or "recover" in o.lower()):
        return ("Missing checkpoint", "")
    m = re.search(r"unused parameter[s]?[^\n]*?'([^']+)'", o, re.I)
    if m: return ("Unused param", m.group(1))
    m = re.search(r"\[([^\]\n]+)\][^\n]*\bunused\b", o)
    if m: return ("Unused param", m.group(1))
    if "were unused" in o or "was unused" in o:
        m = re.search(r"([A-Za-z0-9_]+/[A-Za-z0-9_/]+)", o)
        return ("Unused param", m.group(1) if m else "")
    m = re.search(r"no variable '([^']+)' found for (?:substitution|use)", o)
    if m:
        return ("Unset CLI variable", m.group(1))
    if "CheckIntegrityAction" in o:
        return ("Integrity action", "")
    if "MooseEnum" in o and ("Invalid" in o or "invalid" in o):
        return ("Invalid value", "")
    if re.search(r"not a subset", o):
        return ("Block subset error", "")
    if '"Executioner" does not exist' in o:
        return ("Missing Executioner", "")
    if re.search(r"Variable '([^']+)'.*does not exist|does not exist.*Variable", o):
        return ("Missing variable", "")
    # Other: capture the first meaningful error line as a short snippet
    m = re.search(r"\*\*\* ERROR \*\*\*\s*\n(.*?)(?:\n|$)", o, re.S)
    snippet = ""
    if m:
        for ln in m.group(1).splitlines():
            ln = ln.strip()
            ln = re.sub(r"^/\S+\.i:[\d.]+:\s*", "", ln)  # strip file:line prefix
            if ln and not ln.startswith("The following"):
                snippet = ln
                break
    if len(snippet) > 70:
        snippet = snippet[:67] + "..."
    return ("Other", snippet)

def check_one(arg):
    app, path = arg
    exe = APPS[app]
    d, fn = os.path.dirname(path), os.path.basename(path)
    try:
        r = subprocess.run([exe, "-i", fn, "--check-input"],
                           cwd=d, capture_output=True, text=True, timeout=120)
        out = r.stdout + r.stderr
        ok = (r.returncode == 0) and ("Syntax OK" in out or "Parser/Factory" in out)
        if r.returncode == 0:
            ok = True
    except subprocess.TimeoutExpired:
        return (app, path, False, ("Other", "timeout"))
    except Exception as e:
        return (app, path, False, ("Other", str(e)[:40]))
    if ok:
        return (app, path, True, ("", ""))
    return (app, path, False, categorize(out))

def tex_esc(s):
    s = re.sub(r"[\x00-\x1f\x7f]", "", s)
    return (s.replace("\\", r"\textbackslash{}").replace("_", r"\_")
             .replace("&", r"\&").replace("%", r"\%").replace("#", r"\#")
             .replace("$", r"\$").replace("{", r"\{").replace("}", r"\}"))

def main():
    files = discover()
    print(f"Discovered {len(files)} input files; running --check-input ...", flush=True)
    with ProcessPoolExecutor(max_workers=8) as ex:
        results = list(ex.map(check_one, files))
    print("done.", flush=True)

    # group by app / first-level project dir
    groups = {}
    for app, path, ok, cat in results:
        rel = os.path.relpath(path, os.path.join(FEECM, "projects", app))
        parts = rel.split(os.sep)
        proj = parts[0]
        sub = os.sep.join(parts[1:])
        groups.setdefault((app, proj), []).append((sub, ok, cat))

    total = len(results)
    npass = sum(1 for r in results if r[2])
    nfail = total - npass

    # per-app
    app_stats = {}
    for app, path, ok, cat in results:
        s = app_stats.setdefault(app, [0, 0])
        s[0] += 1
        if ok: s[1] += 1

    # failure categories
    catcount = {}
    classcount = {}
    for app, path, ok, cat in results:
        if ok: continue
        catcount[cat[0]] = catcount.get(cat[0], 0) + 1
        if cat[0] == "Missing class" and cat[1]:
            classcount[cat[1]] = classcount.get(cat[1], 0) + 1

    out = []
    W = out.append
    today = date.today().isoformat()
    W(r"\documentclass[10pt,letterpaper]{article}")
    W(r"\usepackage[margin=1in]{geometry}")
    W(r"\usepackage{longtable}")
    W(r"\usepackage{booktabs}")
    W(r"\usepackage{xcolor}")
    W(r"\usepackage{hyperref}")
    W(r"\usepackage{array}")
    W(r"\hypersetup{colorlinks=true, linkcolor=blue, urlcolor=blue}")
    W(r"\setlength{\tabcolsep}{4pt}")
    W(r"\renewcommand{\arraystretch}{1.1}")
    W(r"\definecolor{passgreen}{RGB}{0,120,0}")
    W(r"\definecolor{failred}{RGB}{180,0,0}")
    W(r"\newcommand{\PASS}{\textcolor{passgreen}{\textbf{PASS}}}")
    W(r"\newcommand{\FAIL}{\textcolor{failred}{\textbf{FAIL}}}")
    W(r"\newcommand{\code}[1]{\texttt{\small #1}}")
    W(r"\title{FEECM Project Audit Report}")
    W(r"\author{Per-project \texttt{-{}-check-input} verification across \texttt{feecm/projects/}}")
    W(r"\date{%s}" % today)
    W(r"\begin{document}")
    W(r"\maketitle")
    W("")
    W(r"\section{Scope and method}")
    W(r"Every input file under \code{feecm/projects/<app>/} was run through the "
      r"corresponding app executable with \code{-{}-check-input}. This validates "
      r"input syntax, parser/factory registrations, mesh/file references, and "
      r"basic integrity, but does \emph{not} time-step the simulation. An input "
      r"that passes \code{-{}-check-input} may still fail at solve time; a "
      r"failing input is definitely broken at parse/setup.")
    W("")
    W(r"As of this audit the repository has \textbf{two} apps (the former "
      r"\code{ec\_beta} and \code{eel} apps were consolidated into \code{ecm\_test} "
      r"and \code{tecm\_test} respectively and removed):")
    W(r"\begin{itemize}")
    W(r"  \item \code{ecmTest} $\to$ \code{ecm\_test/ecm-opt}")
    W(r"  \item \code{tecmTest} $\to$ \code{tecm\_test/tecm\_test-opt} "
      r"(now also hosts the former \code{projects/eel/} work)")
    W(r"\end{itemize}")
    W("")
    W(r"\section{Headline results}")
    W(r"\textbf{Total inputs audited:} %d\\" % total)
    W(r"\textbf{PASS:} %d (\textcolor{passgreen}{%.1f\%%})\\" % (npass, 100*npass/total))
    W(r"\textbf{FAIL:} %d (\textcolor{failred}{%.1f\%%})\\" % (nfail, 100*nfail/total))
    W("")
    W(r"\subsection{Per-app summary}")
    W(r"\begin{tabular}{lrrrr}")
    W(r"\toprule")
    W(r"App & Inputs & PASS & FAIL & Pass rate \\")
    W(r"\midrule")
    for app in ("ecmTest", "tecmTest"):
        t, p = app_stats.get(app, [0, 0])
        rate = 100*p/t if t else 0
        W(r"\code{%s} & %d & %d & %d & %.1f\%% \\" % (app, t, p, t-p, rate))
    W(r"\bottomrule")
    W(r"\end{tabular}")
    W("")
    W(r"\section{Failure categories}")
    W(r"The %d failures break down as follows:" % nfail)
    W(r"\begin{longtable}{lr}")
    W(r"\toprule")
    W(r"Category & Count \\")
    W(r"\midrule")
    W(r"\endhead")
    for c, n in sorted(catcount.items(), key=lambda x: -x[1]):
        W(r"\code{%s} & %d \\" % (tex_esc(c), n))
    W(r"\bottomrule")
    W(r"\end{longtable}")
    if classcount:
        W(r"\subsection{Most-requested missing classes}")
        W(r"\begin{tabular}{lr}")
        W(r"\toprule")
        W(r"Class name & Inputs blocked \\")
        W(r"\midrule")
        for c, n in sorted(classcount.items(), key=lambda x: -x[1]):
            W(r"\code{%s} & %d \\" % (tex_esc(c), n))
        W(r"\bottomrule")
        W(r"\end{tabular}")
    W("")
    W(r"\section{Per-project results}")
    W(r"For each project the table lists every audited input with its status "
      r"and (on failure) a short categorized issue.")
    W("")

    cur_app = None
    for (app, proj) in sorted(groups):
        if app != cur_app:
            cur_app = app
            W(r"\subsection{App: \code{%s}}" % tex_esc(app))
        rows = sorted(groups[(app, proj)])
        gp = sum(1 for _, ok, _ in rows if ok)
        W(r"\subsubsection{\code{%s/%s}\hfill (%d/%d pass)}" %
          (tex_esc(app), tex_esc(proj), gp, len(rows)))
        if gp == len(rows):
            names = ", ".join(r"\code{%s}" % tex_esc(s) for s, _, _ in rows)
            W(r"All %d input(s) pass. Files: %s." % (len(rows), names))
            W("")
            continue
        W(r"\begin{longtable}{p{0.46\textwidth} l p{0.36\textwidth}}")
        W(r"\toprule")
        W(r"Input & Status & Issue \\")
        W(r"\midrule")
        W(r"\endhead")
        for sub, ok, cat in rows:
            if ok:
                W(r"\code{%s} & \PASS & --- \\" % tex_esc(sub))
            else:
                issue = tex_esc(cat[0])
                if cat[1]:
                    issue += r": \code{%s}" % tex_esc(cat[1])
                W(r"\code{%s} & \FAIL & %s \\" % (tex_esc(sub), issue))
        W(r"\bottomrule")
        W(r"\end{longtable}")
        W("")

    W(r"\end{document}")

    texpath = os.path.join(FEECM, "projects", "doc", "project_audit.tex")
    with open(texpath, "w") as f:
        f.write("\n".join(out) + "\n")
    print(f"Wrote {texpath}")
    print(f"TOTAL {total}  PASS {npass}  FAIL {nfail}")
    for app in ("ecmTest", "tecmTest"):
        t, p = app_stats.get(app, [0, 0])
        print(f"  {app}: {p}/{t} pass")

if __name__ == "__main__":
    main()
