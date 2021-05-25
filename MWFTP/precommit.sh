# Checks to run on the code before commiting.

set -e

cd $(dirname $0)

# py_files includes only non-test files.
py_files="$(find | grep -P '\.py$' | grep -vP 'test_.*\.py$')"
test_files="$(find | grep -P 'test_.*\.py$')"

# Run all tests and print coverage.
echo -e '\e[1;34mTesting\e[0m'
python -m coverage run -m unittest
python -m coverage report --include "./*" --omit "test_*"

# Lint files.
# We use different parameters for tests and for regular scripts (for example, we
# don't require comments on tests).
#
# To view details of disabled checks, run:
# pylint --help-msg $disabled_checks_tests
disabled_checks=C0103,W0511,R0913
disabled_checks_tests=$disabled_checks,C0114,C0115,C0116,R0201

echo -e '\e[1;34mLinting\e[0m'
pylint -d $disabled_checks $py_files
pylint -d $disabled_checks_tests $test_files
flake8 --max-line-length=100

echo -e '\e[1;34mSuccess!\e[0m'
