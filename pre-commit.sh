#!/usr/bin/env bash

black ./
RESULT=$?
if [[ ${RESULT} -ne 0 ]]; then
	printf "Failed to format python scripts with black. Commit rejected"
	exit ${RESULT}
fi

shfmt -l -w ./
RESULT=$?

if [[ ${RESULT} -ne 0 ]]; then
	printf "Failed to format bash scripts using shfmt. Commit rejected"
fi
exit ${RESULT}
