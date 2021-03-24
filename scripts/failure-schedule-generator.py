#!/usr/bin/env python3

import math

RANKS_PER_NODE = 20
MAX_FAILED_NODES = 10
FAIL_EVERY = 1000000

for step in range(1, MAX_FAILED_NODES + 1):
    with open(f'failure_{step:02}.schedule', 'w') as schedule:
        call_no = FAIL_EVERY
        failed_nodes = 0
        while failed_nodes <= MAX_FAILED_NODES:
            for node in range(failed_nodes, failed_nodes + step):
                for rank in range(node * RANKS_PER_NODE, (node + 1) * RANKS_PER_NODE):
                    schedule.write(f'{rank} {call_no}\n')
                failed_nodes += 1
            call_no += FAIL_EVERY


