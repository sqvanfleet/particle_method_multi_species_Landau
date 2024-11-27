#!/bin/bash

nice -n14 ionice -c2 -n7 matlab -batch "Example_6_n_50_dt_01_975"
process_id=$!
echo "PID: $process_id"
wait $process_id
echo "Exit status: $?"
nice -n14 ionice -c2 -n7 matlab -batch "Example_6_n_50_dt_01_13"
process_id=$!
echo "PID: $process_id"
wait $process_id
echo "Exit status: $?"

