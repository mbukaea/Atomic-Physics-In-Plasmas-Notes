#!/usr/bin/env bash

./run_tests.sh

RESULT=$?

if [[ ${RESULT} -ne 0 ]]
then
    RED='\033[0;31m'
    printf  "${RED}Unit tests failed. Commit rejected"
fi
exit ${RESULT}
