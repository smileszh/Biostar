#!/bin/bash
#
# Report timing information for a command.
#
# Usage:
#
#    time.sh command
#
# Generates both an output on the screen and saves the timing information into a log file see $LOG variable.
#
# Uses the GNU gtime command
#
# On Linux create a gtime link to time with
#
#    ln -s /usr/bin/time ~/bin/gtime
#
# On Mac OS you get gtime with:
#
#   brew install gnu-time
#

# The name of the logfile.
LOG=timelog.txt

# The name of the timing command to use.
TIME=gtime

# Check if the timer is installed.
if ! command -v ${TIME} &> /dev/null; then
  echo "#"
  echo "# Timing command not found: ${TIME}"
  echo "#"
  echo "# Linux setup: ln -s /usr/bin/time ~/bin/gtime"
  echo "#"
  echo "# MacOS setup: brew install gnu-time"
  echo "#"
  exit 1
fi

# Check if the command is passed as parameter
if [ -z "$1" ]; then
  echo "#"
  echo "# Usage: $0 <command to time>"
  echo "#"
  exit 1
fi

# Bash strict mode
set -ue

# Print the command that is being timed to the stderr.
>&2 echo "# running: $@"

# Run the command line passed to this script.
# %e   elapsed real time (wall clock) in seconds
# %E   elapsed real time (wall clock) in [hour:]min:sec
# %M   maximum resident set size in KB
# %x   exit status of the command
# %C   name and command line arguments of the command being timed
# -a   append the output to the file
# -o   save the output to the file

# Format the output as a JSON object.
DATE=$(date '+%Y-%m-%d, %H:%M:%S')

# Store the hostname
HOST=$(hostname)
${TIME} -q -f "{ \"code\":\"%x\", \"sec\":\"%e\", \"time\":\"%E\", \"mem\":\"%M\", \"cmd\":\"%C\", \"date\":\"$DATE\", \"host\":\"$HOST\" }" -a -o ${LOG} $@

# If the LOG file exists
if [ -f "$LOG" ]; then
  # Print the last line of the log file to the standard error.
  echo "# ${LOG}: $(tail -1 ${LOG})" >&2
fi

