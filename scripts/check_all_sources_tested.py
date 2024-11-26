#!/usr/bin/env python3

from glob import glob
from generate_include_graph import construct_edges

if __name__ == "__main__":
    # We first construct a list of all .hpp files in src/
    sources = ["src/"+name for name in glob('**/*.hpp', root_dir="src", recursive=True)]
    # and a list of all include statements in src/ and testsuite/
    includes = construct_edges(base_dir=".", quiet=True)

    for include in includes:
        try:
            # We then remove all #included files from the list
            sources.remove(include[1])
        except ValueError:
            continue

    # and print all remaining ones
    print("The following source-files are not being tested anywhere")
    for s in sources:
        print(s)

    #TODO: We should - once all tests are implemented - assert tht the list is empty.
