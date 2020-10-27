#!/bin/bash
pbossid=$(ps uax|grep -m1 -i "/Applications/Password Boss.app/Contents/MacOS/Password Boss"|grep -v grep|awk '{print $2}')
if [ "$pbossid" != "" ];then
kill -9 $pbossid
fi
