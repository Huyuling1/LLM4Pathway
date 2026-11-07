"""Render English Beamer tables/figures directly from deterministic audit outputs."""

import argparse
import csv
import gzip
import json
import math
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from reactome_pipeline.common import file_sha, read_json, read_jsonl, write_bytes, write_json, write_jsonl


def tex(value):
    substitutions = {"\\": r"\textbackslash{}", "&": r"\&", "%": r"\%", "$": r"\$", "#": r"\#",
                     "_": r"\_", "{": r"\{", "}": r"\}", "~": r"\textasciitilde{}", "^": r"\textasciicircum{}"}
    return "".join(substitutions.get(c, c) for c in str(value))


def number(n):
    if isinstance(n, bool):
        return "Yes" if n else "No"
    return f"{n:,}" if isinstance(n, int) else f"{n:.2f}" if isinstance(n, float) else str(n)


def table(headers, rows, widths=None):
    columns = widths or "l" + "r" * (len(headers) - 1)
    lines = [r"\begin{tabularx}{\textwidth}{" + columns + "}", r"\toprule",
             " & ".join(map(tex, headers)) + r" \\", r"\midrule"]
    lines += [" & ".join(tex(number(x)) for x in row) + r" \\" for row in rows]
    return "\n".join(lines + [r"\bottomrule", r"\end{tabularx}"])


def frame(title, body, note=None):
    return (r"\begin{frame}{" + tex(title) + "}\n" + body +
            ("\n" + r"\auditnote{" + tex(note) + "}" if note else "") + "\n" + r"\end{frame}" + "\n")


def plot(histograms, xlabel, ylabel, width="12cm", height="4.3cm"):
    lines = [r"\begin{tikzpicture}\begin{axis}[ybar stacked,bar width=15pt,width=" + width + ",height=" + height +
             r",ymin=0,enlarge x limits=0.12,xtick=data,axis lines*=left,tick label style={font=\scriptsize},label style={font=\small},"
             + "xlabel={" + tex(xlabel) + "},ylabel={" + tex(ylabel) + "},legend style={font=\\scriptsize,at={(0.5,1.12)},anchor=south,legend columns=-1},ytick distance=1]"]
    keys = sorted({int(k) for _, h in histograms for k in h})
    for (label, hist), colour in zip(histograms, ["accent", "warm"]):
        coordinates = " ".join(f"({key},{hist.get(str(key), 0)})" for key in keys)
        lines += [r"\addplot[fill=" + colour + ",draw=" + colour + "] coordinates {" + coordinates + "};",
                  r"\addlegendentry{" + tex(label) + "}"]
    return "\n".join(lines + [r"\end{axis}\end{tikzpicture}"])


def export_csv(path, headers, rows):
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(headers)
        writer.writerows(rows)


def build(cache_report, case, out, generated, extraction_suite_report=None):
    case, out, generated = Path(case), Path(out), Path(generated)
    out.mkdir(parents=True, exist_ok=True)
    generated.mkdir(parents=True, exist_ok=True)
    cache = read_json(cache_report)
    report = read_json(case / "audit/report.json")
    review = read_json(case / "review/summary.json")
    manifest = read_json(case / "snapshot.json")
    reactions = {r["id"]: r for r in read_jsonl(case / "normalized/reactions.jsonl")}
    edges = list(read_jsonl(case / "graph/edges.jsonl"))
    chains = [r for name in ("complete_chains.jsonl", "cycle_stopped_chains.jsonl") for r in read_jsonl(case / "paths" / name)]
    aliases = {rid: "R" + str(i) for i, rid in enumerate(sorted(reactions), 1)}
    if report["graph"]["graph_sha256"] != report["paths"]["graph_sha256"]:
        raise ValueError("Case graph and chains do not match")
    frames = []
    scope_rows = [
        ("Expected reactions", cache["normalization"]["expected"], manifest["expected_reactions"]),
        ("Available reactions", cache["graph"]["reactions"], report["graph"]["reactions"]),
        ("Evidence edges", cache["graph"]["evidence_edges"], report["graph"]["evidence_edges"]),
        ("Reaction pairs", cache["graph"]["reaction_pairs"], report["graph"]["reaction_pairs"]),
        ("Input complete", str(cache["graph"]["input_complete"]), str(report["graph"]["input_complete"])),
        ("Enumeration complete", str(cache["paths"]["enumeration_complete"]), str(report["paths"]["enumeration_complete"]))]
    frames.append(frame("Two audit populations", table(["Measure", "Cached global subset", "Complete pathway case"], scope_rows, "Xrr"),
        "Case completeness refers to the pathway membership list. Its internal graph is an audit view; global matching remains unrestricted."))
    frames.append(frame("Cached graph: composition and topology", r"\small" + table(
        ["Measure", "Count"], [("Species: " + name, count) for name, count in cache["species"].items()] +
        [("Branching reactions", cache["topology"]["branching_nodes"]), ("Merging reactions", cache["topology"]["merging_nodes"]),
         ("Isolated reactions", cache["topology"]["isolated_nodes"]), ("Cyclic strongly connected components", cache["topology"]["cyclic_components"]),
         ("Largest cyclic component (reactions)", cache["topology"]["largest_cyclic_component"])], "Xr"),
        "The available cache follows inventory download order. These counts describe this cache, not the whole Reactome database."))
    frames.append(frame("Cached graph: typed evidence", table(["Type", "Evidence edges"], sorted(cache["edge_types"].items()), "Xr"),
        "One evidence edge is (source reaction, target reaction, role type, physical-entity ID). Multiple evidence edges can share a reaction pair."))
    p = cache["paths"]
    frames.append(frame("Large-graph enumeration: observed checkpoint", table(["Measure", "Observed"], [
        ("DFS operations", p["dfs_steps"]), ("Candidate terminal paths", p["candidate_leaves"]),
        ("Emitted sink chains", p["sink_chains"]), ("Emitted cycle-stopped chains", p["cycle_stopped_chains"]),
        ("Eligible starts entered", str(p["starts_entered"]) + " / " + str(p["starts"])),
        ("Enumeration complete", str(p["enumeration_complete"]))], "Xr"),
        "Run paused with its checkpoint saved. A terminal path is rejected if its first reaction still has an unvisited predecessor. This DFS prefix is not a representative chain sample."))
    frames.append(frame("Complete case: glutathione synthesis and recycling",
        tex(manifest["pathway_id"]) + r"\par\vspace{0.3cm}" + table(["Measure", "Result"], [
        ("Member reactions downloaded", str(review["present_reactions"]) + " / " + str(review["expected_reactions"])),
        ("Species", ", ".join(report["species"])), ("Inferred reactions", report["inference"].get("True", 0)),
        ("Internal evidence edges", len(edges)), ("Reactions in cycles", report["topology"]["reactions_in_cycles"]),
        ("Boundary evidence in cached context", review["context_boundary_evidence_edges"])], "Xr"),
        "Membership is frozen from ContentService containedEvents. The case contains branching, recycling, negative feedback and distinct compartments."))

    active = [rid for rid in reactions if any(e["source"] == rid or e["target"] == rid for e in edges)]
    isolated = [rid for rid in reactions if rid not in active]
    diagram = [r"\begin{columns}[T]\begin{column}{0.53\textwidth}",
               r"\begin{tikzpicture}[>=Stealth,every node/.style={font=\scriptsize},scale=0.9]"]
    for i, rid in enumerate(active):
        angle = math.radians(90 - i * 360 / len(active))
        x, y = 2.35 * math.cos(angle), 2.05 * math.sin(angle)
        diagram.append(f"\\node[draw=accent,fill=white,circle,minimum size=0.65cm] ({aliases[rid]}) at ({x:.3f},{y:.3f}) {{{aliases[rid]}}};")
    for i, rid in enumerate(isolated):
        diagram.append(f"\\node[draw=muted,fill=white,circle,minimum size=0.65cm] ({aliases[rid]}) at (3.1,{-1.8+i:.2f}) {{{aliases[rid]}}};")
    for edge in edges:
        a, b = aliases[edge["source"]], aliases[edge["target"]]
        colour = "warm" if edge["kind"] == "negative_regulation" else "accent"
        style = "loop right,looseness=6" if a == b else "bend left=10"
        diagram.append(f"\\draw[->,{colour},thick] ({a}) to[{style}] ({b});")
    diagram += [r"\end{tikzpicture}\end{column}\begin{column}{0.45\textwidth}\scriptsize",
                table(["Label", "Reaction name"], [(aliases[rid], reactions[rid]["display_name"]) for rid in sorted(reactions)], "lX"),
                r"\end{column}\end{columns}"]
    frames.append(frame("Exact-ID graph of the complete case", "\n".join(diagram),
        "Blue: material flow. Orange: negative regulation. R7 uses extracellular entities; its same-name cytosolic counterparts have different stable IDs."))

    for start in range(0, len(reactions), 3):
        ids = sorted(reactions)[start:start + 3]
        frames.append(frame("Recorded reaction descriptions" + (" (continued)" if start else ""),
            r"\footnotesize" + table(["Label", "Fixed-template text"], [(aliases[rid], reactions[rid]["text"]) for rid in ids], "lX"),
            "Names, compartments, coefficients and catalyst/regulator labels are taken from recorded fields."))

    metric = report["chain_metrics"]
    frames.append(frame("Complete-case chains: termination and length", table(["Termination", "Chains", "Fraction"], [
        ("Sink-ending", metric["sink_chains"], f"{100 * metric['sink_fraction']:.0f}%"),
        ("Cycle-stopped", metric["cycle_stopped_chains"], f"{100 * metric['cycle_stopped_fraction']:.0f}%")], "Xrr") +
        "\n" + plot([(label, metric["lengths"].get(key, {}).get("histogram", {})) for label, key in [("Sink", "sink"), ("Cycle-stopped", "cycle_blocked")]], "Reactions per chain", "Chains", height="3.4cm"),
        f"All {metric['total_chains']} maximal simple chains enumerated. Mean length {metric['lengths']['all']['mean']:.1f}; median {metric['lengths']['all']['median_nearest_rank']}; range {metric['lengths']['all']['min']}–{metric['lengths']['all']['max']}. The only sink chain is isolated R7."))
    coverage = {r["reaction_id"]: r["chain_occurrences"] for r in read_jsonl(case / "audit/reaction_chain_coverage.jsonl")}
    c = metric["reaction_coverage"]
    frames.append(frame("Reaction reuse across complete-case chains", table(["Reaction", "Chain appearances"], [(aliases[rid] + ": " + rid, coverage[rid]) for rid in sorted(reactions)], "Xr"),
        f"{c['total_appearances']} appearances across {c['covered_reactions']} distinct reactions; {c['reactions_in_multiple_chains']} appear in multiple chains. Repeat-appearance fraction = (appearances − distinct reactions) / appearances = {100*c['repeat_appearance_fraction']:.1f}%."))

    # Stable choice: longest chain, then lexical reaction-ID sequence.
    chosen = sorted(chains, key=lambda row: (-len(row["reaction_ids"]), row["reaction_ids"]))[0]
    write_json(out / "example_chain.json", chosen)
    def matched_label(item):
        sources = list(dict.fromkeys(aliases[p["reaction_id"]] for p in item["providers"]))
        return item["display_name"] + " (" + ", ".join(sources) + ")"
    prefix_rows = [(aliases[step["reaction_id"]], "; ".join(matched_label(x) for x in step["covered_inputs"]) or "None",
                    "; ".join(x["display_name"] for x in step["uncovered_inputs"]) or "None") for step in chosen["input_coverage"]]
    frames.append(frame("Worked chain: prefix input coverage", tex(" → ".join(aliases[rid] for rid in chosen["reaction_ids"])) +
        r"\par\vspace{0.15cm}\scriptsize" + table(["Step", "Matched inputs (earlier source)", "Uncovered inputs"], prefix_rows, "lXX"),
        "The displayed chain is selected by longest length, then stable-ID order. Earlier sources may be non-adjacent. ATP and water remain visible."))
    closing = [(aliases[e["source"]] + " → " + aliases[e["target"]],
                "; ".join(x["kind"] + ": " + x["entity_name"] for x in e["evidence"])) for e in chosen["closing_edges"]]
    frames.append(frame("Worked chain: why it stops", table(["Closing connection", "Preserved evidence"], closing, "lX") +
        r"\vspace{0.4cm}\par " + tex("All outgoing targets of the last reaction already occur in the chain. The record is cycle-stopped."),
        "The complete graph retains these edges. Repeating a reaction is disallowed only during chain enumeration."))

    cur = cache["currency"]
    frames.append(frame("Currency exclusions: effect on the cached graph", table(["Measure", "Count"], [
        ("Excluded physical entities", cur["excluded_entities"]), ("Excluded evidence edges", cur["excluded_evidence_edges"]),
        ("Reactions losing all upstream connections", len(cur["reactions_losing_all_upstream"])),
        ("Reactions losing all downstream connections", len(cur["reactions_losing_all_downstream"])),
        ("Reactions becoming isolated", len(cur["reactions_becoming_isolated"]))], "Xr"),
        "Before/after comparison uses the same cached reactions. Endpoint losses follow removal of the entire frozen currency list. Self-loops count as connections."))
    currencies = sorted(cur["molecules"], key=lambda x: (-x["excluded_evidence_edges"], x["identifier"]))
    currency_rows = [(x["name"], x["identifier"], x["excluded_evidence_edges"], len(x["affected_reactions_losing_all_upstream"]),
                      len(x["affected_reactions_losing_all_downstream"]), len(x["affected_reactions_becoming_isolated"])) for x in currencies]
    export_csv(out / "currency_molecules.csv", ["molecule", "identifier", "excluded_evidence", "affected_lost_upstream", "affected_lost_downstream", "affected_newly_isolated"], currency_rows)
    for offset in range(0, len(currency_rows), 7):
        frames.append(frame(f"Currency audit by identifier ({offset//7 + 1}/{math.ceil(len(currency_rows)/7)})",
            r"\scriptsize" + table(["Molecule", "ChEBI ID", "Edges", "Lost in", "Lost out", "Isolated"], currency_rows[offset:offset + 7], "Xlrrrr"),
            "Lost in/out/isolated: reactions touched by this molecule among all-list endpoint losses. These sets can overlap across molecules; do not sum their counts."))
    example_rows = []
    for molecule in currencies[:2]:
        if molecule["examples"]:
            e = molecule["examples"][0]
            example_rows.append((molecule["name"], e["source"] + " → " + e["target"], e["entity_name"], str(e["pair_retained_by_other_evidence"])))
    frames.append(frame("Inspecting excluded connections", r"\scriptsize" + table(
        ["Molecule", "Reaction pair", "Shared entity", "Pair retained?"], example_rows, "lXXl") +
        r"\par\vspace{0.45cm}\normalsize " + tex(f"In the complete glutathione case, {report['currency']['excluded_entities']} currency entities occur, but they form {report['currency']['excluded_evidence_edges']} candidate evidence edges. The frozen exclusion rule removes no internal connection in this case."),
        "Pair retained means another non-currency evidence edge still connects the same reactions. All per-entity examples and endpoint IDs are saved in the audit JSON."))

    preceding_rows = [(label, cache["preceding_event"][key], report["preceding_event"][key]) for label, key in [
        ("Recorded precedingEvent pairs", "pairs_within_graph"),
        ("Shared by both methods", "overlap"),
        ("Entity matching only", "only_entity_matching"),
        ("precedingEvent only", "only_preceding"),
        ("Matched negativePrecedingEvent records", "matched_negative_preceding")]]
    export_csv(out / "preceding_event_comparison.csv", ["metric", "cached_graph", "complete_case"], preceding_rows)
    frames.append(frame("Audit against recorded preceding events",
        r"\small" + table(["Measure", "Cached graph", "Complete case"], preceding_rows, "Xrr"),
        "Compare directed reaction pairs with both endpoints inside each graph. Multiple evidence types collapse to one pair. negativePrecedingEvent matches count explicit rejection records; external events are reported separately."))

    preceding_pairs = {(event["id"], rid) for rid, reaction in reactions.items()
                       for event in reaction["preceding_events"] if event["id"] in reactions}
    matched_pairs = {(edge["source"], edge["target"]) for edge in edges}
    if (len(matched_pairs & preceding_pairs), len(matched_pairs - preceding_pairs), len(preceding_pairs - matched_pairs)) != (
            report["preceding_event"]["overlap"], report["preceding_event"]["only_entity_matching"], report["preceding_event"]["only_preceding"]):
        raise ValueError("Case precedingEvent audit does not match its reaction and edge records")
    extra_edges = [edge for edge in edges if (edge["source"], edge["target"]) not in preceding_pairs]
    write_jsonl(out / "case_entity_only_edges.jsonl", extra_edges)
    extra_rows = [(aliases[edge["source"]] + " → " + aliases[edge["target"]], edge["kind"], edge["entity_name"])
                  for edge in extra_edges]
    supported = sum((step["source"], step["target"]) in preceding_pairs for step in chosen["transitions"])
    frames.append(frame("Complete case: connections beyond precedingEvent",
        r"\small" + table(["Connection", "Matching rule", "Shared entity"], extra_rows, "lXX") +
        r"\par\vspace{0.4cm}" + tex("Displayed chain: " + " → ".join(aliases[rid] for rid in chosen["reaction_ids"])) +
        r"\par\vspace{0.2cm}" + tex(f"{supported} / {len(chosen['transitions'])} forward connections also have a precedingEvent entry."),
        "An absent precedingEvent entry is an unmatched candidate, not an explicit rejection. The graph retains both additional connections for review."))

    if "extraction_audit" in cache and "extraction_audit" in report:
        source_audit, case_source = cache["extraction_audit"], report["extraction_audit"]
        if not cache.get("extraction_matches_current_normalized") or not report.get("extraction_matches_current_normalized"):
            raise ValueError("Extraction audit is stale relative to normalized data")
        frames.append(frame("Extraction audit: raw slots versus normalized roles", table(
            ["Measure", "Cached reactions", "Complete case"],
            [(label, source_audit[key], case_source[key]) for label, key in [
                ("Reactions compared", "reactions_compared"), ("Role groups checked", "role_checks"),
                ("Role groups matching source", "matched_role_checks"),
                ("Participants in raw role slots", "source_participants"),
                ("Participants in normalized roles", "normalized_participants"),
                ("Source comparison passed", "source_comparison_passed")]], "Xrr"),
            "The auditor traverses raw root-role slots without calling the normalizer. Checks cover counts, pointers, entity IDs, names, classes, coefficients, compartments, ChEBI references and catalyst/regulation fields."))
        empty_rows = [(role, sum(reasons.values()), "; ".join(f"{reason}: {count}" for reason, count in sorted(reasons.items())) or "None")
                      for role, reasons in sorted(source_audit["empty_roles"].items())]
        frames.append(frame("Empty roles: trace each absence to the source", r"\small" + table(
            ["Normalized role", "Empty reactions", "Observed source state"], empty_rows, "l r X"),
            "These empty roles have no corresponding field in the saved response. The audit distinguishes absent, null, empty-list, actor-free group, unresolved reference and extraction omission states. It does not infer a biological reason for an absent field."))
        write_json(out / "extraction-cache-audit.json", source_audit)
        write_json(out / "extraction-case-audit.json", case_source)
        if extraction_suite_report:
            suite = read_json(extraction_suite_report)
            suite_classes = {row["schema_class"]: row["compared"] for row in suite["schema_class_coverage"]}
            class_rows = [(row["schema_class"], row["inventory"], row["compared"], suite_classes.get(row["schema_class"], 0))
                          for row in source_audit["schema_class_coverage"]]
            frames.append(frame("Real-record coverage across event subclasses", r"\small" + table(
                ["Event subclass", "Full inventory", "Cache checked", "Frozen examples"], class_rows, "Xrrr") +
                r"\par\vspace{0.3cm}" + tex(f"Example suite: {suite['reactions_compared']} reactions, {suite['source_participants']} participants; source comparison passed: {suite['source_comparison_passed']}."),
                "One record per inventory subclass: cached records first, then stable-ID order. Missing examples are fetched with a release check. Selection, raw responses and checksums are frozen for offline reruns; coverage is one example per class."))
            write_json(out / "extraction-suite-audit.json", suite)

    frames.append(frame("Source checks and human review", table(["Check", "Result"], [
        ("Pathway reaction membership", str(review["present_reactions"]) + " / " + str(review["expected_reactions"])),
        ("Participant pointers matching source stable IDs", str(review["source_pointer_passes"]) + " / " + str(review["source_pointer_checks"])),
        ("Reaction review items", len(reactions)), ("Edge review items", len(edges)),
        ("Human review status", "; ".join(f"{k}: {v}" for k, v in sorted(review["manual_review_status_counts"].items())))], "Xr"),
        "Pointer checks verify extracted participant identity. Human review covers field completeness, stoichiometry, edge roles and biological interpretation; record decisions in manual_review.csv."))
    frames.append(frame("Files for reproducing and extending this audit", r"\small" + table(["Artifact", "Purpose"], [
        ("examples/pathway_snapshot/", "Complete raw case, membership responses and checksums"),
        ("reactome_pipeline/coverage.py", "All-prefix exact input matching"),
        ("reactome_pipeline/audit_details.py", "Chain and currency statistics"),
        ("reactome_pipeline/review.py", "Pathway capture and source-linked review packet"),
        ("reports/demo-audit/", "Frozen audit records, tables and example chain"),
        ("slides/main.tex", "English Beamer source"),
        ("scripts/build_audit_slides.py", "Deterministic tables and figures")], "lX"),
        "Rebuild commands are in slides/README.md. Source hashes and graph/configuration fingerprints are saved with the report."))
    frames.append(frame("Primary source endpoints", r"\footnotesize" + "\n".join([
        r"\url{https://reactome.org/ContentService/data/database/version}\par\vspace{0.3cm}",
        r"\url{https://reactome.org/ContentService/data/schema/ReactionLikeEvent/min}\par\vspace{0.3cm}",
        r"\url{" + manifest["membership_sources"][0]["url"] + r"}\par\vspace{0.3cm}",
        r"\url{" + manifest["membership_sources"][1]["url"] + r"}\par\vspace{0.3cm}",
        r"\texttt{/data/query/enhanced/v2/\{dbId\}}: nested JSON, references and incoming relationships."]),
        "Raw responses and their URLs are preserved in the frozen case. All reaction and entity labels displayed here come from those records."))
    write_bytes(generated / "results.tex", ("% !TeX root = ../main.tex\n" + "\n".join(frames)).encode())
    # A compact frozen cache report supports offline slide regeneration.
    compact = dict(cache)
    compact["normalization"] = {k: v for k, v in cache["normalization"].items() if k != "missing_db_ids"}
    write_bytes(out / "cache-audit.json.gz", gzip.compress((json.dumps(compact, ensure_ascii=False, sort_keys=True, indent=2) + "\n").encode(), mtime=0))
    write_json(out / "case-audit.json", report)
    write_json(out / "review-summary.json", review)
    write_jsonl(out / "case-chains.jsonl", chains)
    write_jsonl(out / "case-reactions.jsonl", reactions.values())
    write_jsonl(out / "case-edges.jsonl", edges)
    write_json(out / "provenance.json", {"reactome_release": manifest["reactome_release"],
        "case_inventory_sha256": manifest["inventory_sha256"], "cache_inventory_sha256": cache["normalization"]["inventory_sha256"],
        "case_graph_sha256": report["graph"]["graph_sha256"], "cache_graph_sha256": cache["graph"]["graph_sha256"],
        "currency_sha256": report["graph"]["currency_sha256"], "chain_record_version": report["paths"]["record_version"],
        "source_hashes": {rid: r["source"]["sha256"] for rid, r in reactions.items()},
        "code_hashes": {p.relative_to(REPO).as_posix(): file_sha(p) for p in sorted((REPO / "reactome_pipeline").glob("*.py"))},
        "report_generator_sha256": file_sha(__file__), "slide_template_sha256": file_sha(REPO / "slides/main.tex")})
    return {"generated_frames": len(frames), "case_chains": len(chains), "output": str(generated / "results.tex")}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache-report", default="data/reactome-98/audit/report.json")
    parser.add_argument("--case", default="data/audit-glutathione-98")
    parser.add_argument("--out", default="reports/demo-audit")
    parser.add_argument("--generated", default="slides/generated")
    parser.add_argument("--extraction-suite-report", default="reports/demo-audit/extraction-suite-audit.json")
    args = parser.parse_args()
    print(json.dumps(build(args.cache_report, args.case, args.out, args.generated,
                           args.extraction_suite_report if Path(args.extraction_suite_report).exists() else None), sort_keys=True))
