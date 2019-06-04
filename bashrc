# get the current directory
CURR_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

export PATH=${CURR_DIR}/platform/bin:$PATH
export LD_LIBRARY_PATH=${CURR_DIR}/platform/lib:$LD_LIBRARY_PATH
