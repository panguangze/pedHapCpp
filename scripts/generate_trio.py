import argparse
import logging
from typing import List
import ped
import ped_utils
from ped_utils import Trio
logger = logging.getLogger(__name__)

def up_to_down(families, out_file):
    for f in families:
        all_trios = ped_utils.get_trios(f)
        top_level_trios = ped_utils.get_top_level_trios(all_trios)
        for t in top_level_trios:
            dad_id = "null"
            mom_id = "null"
            if t.dad != None:
                dad_id = t.dad.id
            if t.mom != None:
                mom_id = t.mom.id
            out_file.write("{} {} {}\n".format(t.child.id, dad_id, mom_id))
        next_level_trios = ped_utils.get_next_level_trios(
            all_trios, top_level_trios)
        while len(next_level_trios) != 0:
            for t in next_level_trios:
                out_file.write("{} {} {}\n".format(t.child.id, t.dad.id, t.mom.id))
            next_level_trios = ped_utils.get_next_level_trios(
                all_trios, next_level_trios)


def down_to_up(families, out_file):
    for f in families:
        all_trios = ped_utils.get_trios(f)
        bottom_level_trios = ped_utils.get_bottom_level_trios(all_trios)
        for t in bottom_level_trios:
            out_file.write("{} {} {}\n".format(t.child.id, t.dad.id, t.mom.id))
        prev_level_trios = ped_utils.get_prev_level_trios(
            all_trios, bottom_level_trios)
        while len(prev_level_trios) != 0:
            for t in prev_level_trios:
                out_file.write("{} {} {}\n".format(t.child.id, t.dad.id, t.mom.id))
            prev_level_trios = ped_utils.get_prev_level_trios(
                all_trios, prev_level_trios)

def main():
    parser = argparse.ArgumentParser("trio phase")
    parser.add_argument(
        '-p', help='pedigree file', required=True, dest='ped_file')
    parser.add_argument(
        '-o', help='phasing order file', required=True, dest='out_file')
    args = parser.parse_args()
    families = ped.open_ped(args.ped_file)
    out_put = open(args.out_file, "w")
    up_to_down(families, out_put)
    out_put.write("# # #\n")
    down_to_up(families, out_put)
if __name__ == "__main__":
    main()
