#!/usr/bin/env python3
# coding: utf-8
#
# Requires paramiko for ssh connections
#
# This script will fetch the list of hostnames this Slurm job is running on, ssh onto them
# and kill the specified MPI processes to simulate core failure for testing.
# We are assuming, that there is one rank per physical node!

import paramiko
import sh
import argparse
import os
import time
import random
import re

# Functions
def parseNodeList(nodelist):
    """Generates a full node list out of the abbreviated version squeue gives"""
    if (re.match("^fh2n[0-9]{4}$", nodelist)):
        return [nodelist]
    else:
        assert(re.match("^fh2n\[([0-9,-]+)\]$", nodelist))
        nodelist = nodelist[5:-1]

        expandedNodeList = []
        for item in nodelist.split(','):
            if '-' in item:
                lower, upper = item.split('-')
                expandedNodeList.extend(range(int(lower), int(upper) + 1))
            else:
                expandedNodeList.append(item)

        return ["fh2n%s" % item for item in expandedNodeList]

# Needed command line tools
import sh
from sh import squeue

# Parse command line arguments
parser = argparse.ArgumentParser(description = 'Simulate MPI rank failure on the ForHLR II cluster.')
parser.add_argument('intervall', type=int, help = 'Intervall in seconds between two failures.')

args = parser.parse_args()
assert(args.intervall > 0)

# Load user credential
with open(os.environ['HOME'] + '/.ssh_pwd') as pwd_file:
    PASSWORD = pwd_file.read()[:-1] # Remove trailing newline
USER = 'co8976'

# Get nodelist
JOB_ID = os.environ['SLURM_JOB_ID']
assert(JOB_ID != '');

nodelist = parseNodeList(str(squeue('--jobs=%s' % JOB_ID, '--format=\%N')).split('\n')[1][1:])
print("[Failure Simulator] Running on nodes " + ' '.join(nodelist))
print("[Failure Simulator] I will send SIGKILL to a random MPI rank every %d second(s)" % args.intervall)

client = paramiko.SSHClient()
client.load_system_host_keys()
client.set_missing_host_key_policy(paramiko.AutoAddPolicy())

while True:
    time.sleep(args.intervall)
    nodeToKill = random.choice(nodelist)
    client.connect(nodeToKill, username=USER, password=PASSWORD)
    stdin, stdout, stderr = client.exec_command('kill -s KILL $(pgrep raxml-ng-mpi | head -1)')
    print("Killed first RAxML rank on node " + nodeToKill)

