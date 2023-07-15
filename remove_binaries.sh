#!/bin/bash
git filter-repo --path-glob '*.jar' --invert-paths --force
git filter-repo --path-glob '*.so' --invert-paths --force
git filter-repo --path-glob '*.pyc' --invert-paths --force
git filter-repo --path-glob '*.pdf' --invert-paths --force
git filter-repo --path-glob '*.png' --invert-paths --force
git filter-repo --path-glob '*.TAB' --invert-paths --force
git filter-repo --path-glob '*.P' --invert-paths --force
git filter-repo --path-glob '*.p' --invert-paths --force
git filter-repo --path-glob '*.DS_Store' --invert-paths --force
git filter-repo --path-glob '*test' --invert-paths --force
git filter-repo --path-glob 'docs/_build/*' --invert-paths --force
git filter-repo --path-glob 'examples/math_examples/*' --invert-paths --force
git filter-repo --path-glob 'examples/numerical_examples/*' --invert-paths --force
