#!/bin/sh
#
# File: append_footer.sh
#
# Usage:
# append_footer.sh [-debug] "name-of-footer-text-file" "input_file" "output_file"
#             usually "~/bin/SparseLab_footer.txt"
# 
# References:
# http://docs.sun.com/source/816-6009-10/channel2.htm#42323
# http://docs.sun.com/source/816-6009-10/channel2.htm#42402
#
# Note: will overwrite Contents.m files, or any file that is exclusively comments
# 
# Modified by Victoria Stodden

if [ "$1" = "-debug" ]
then
shift
set -x
fi

FOOTER_FILE=$1
#FOOTER_FILE=~/bin/ ${FOOTER_FILE}

touch $3

INPUT_FILE=$2
OUTPUT_FILE=$3

TAG="% Standard SparseLab Footer Appended `date`"

cp $INPUT_FILE $OUTPUT_FILE # copy original message part to output destination.

# See if the message was already tagged.
grep "Comments: Standard SparseLab Footer Appended" $MESSAGE_HEADERS >/dev/null
if [ $? -ne 0 ]
then
# add a blank line
echo "" >> $OUTPUT_FILE

# append the disclaimer
cat $FOOTER_FILE >> $OUTPUT_FILE

# Set a directive so the message will be tagged
echo "OUTPUT_DIAGNOSTIC=\"${TAG}\"" > $OUTPUT_OPTIONS
fi

# end script.

