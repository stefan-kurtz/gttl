#!/usr/bin/env python3
"""Generate a graph of all includes in the project using Graphviz"""

import subprocess
from glob import glob


def construct_nodes():
    """Get a list of all C++ Modules in the project"""
    source_files = glob("**/*.[ch]pp", recursive=True)
    return source_files


def construct_edges(base_dir="..", quiet=False):
    """Get a list of all project-internal #include statements"""
    internal_include_lines = (
        subprocess.run(
            ["grep", "-RiIs", "--exclude-dir=.git", "-E", '#include\\W+"', base_dir],
            stdout=subprocess.PIPE,
            check=True,
        )
        .stdout.decode("utf-8")
        .strip()
        .split("\n")
    )

    includes = []

    for line in internal_include_lines:
        file, depends_on = line.split(":")
        depends_on = depends_on.split('"')[1]
        if depends_on.startswith(
            ("alignment", "indexes", "sequences", "threading", "utilities")
        ):
            depends_on = "./src/" + depends_on
        if "/" not in depends_on:
            if not quiet:
                print(
                    "WARNING: NO INCLUDE PATH SPECIFIED IN: "
                    + file
                    + " -> "
                    + depends_on
                )
                print("Resolved by changing to: " + depends_on)
            depends_on = "/".join(file.split("/")[:-1]) + "/" + depends_on
        includes.append((file, depends_on))
    return includes


def main():
    """Primary entry point"""
    # we import this here so that check_all_sources_tested.py does not depend on graphviz.
    # This allows the CI to run it without requiring a graphviz installation
    import graphviz

    # print(json.dumps(DEPENDS, sort_keys=True, indent=4))
    nodes = construct_nodes()
    edges = construct_edges()

    dot = graphviz.Digraph(
        "GTTL",
        comment="Include Relations",
        format="svg",
        engine="neato",
        graph_attr={"overlap": "false", "splines": "true"},
    )

    for node in nodes:
        if node.startswith("tools/"):
            continue
        dot.node(node, node.split("/")[-1])
    for edge in edges:
        if edge[0].startswith("tools/"):
            continue
        if edge[0].startswith("testsuite/"):
            continue
            dot.edge(*edge, constraint="false")
        else:
            dot.edge(*edge)

    dot.render(directory=".", view=False)


if __name__ == "__main__":
    main()
