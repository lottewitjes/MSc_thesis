#! /usr/bin/evn/python

from sys import argv

def select_relevant_hits(input, output):
    hmm_list = ["Chal_sti_synt_C", "Chal_sti_synt_N", "SQHop_cyclase_C", "SQHop_cyclase_N",
                "Terpene_synth_C", "Terpene_synth", "Bet_v_1", "Str_synth", "Dirigent"]
    with open(input, "r") as f1, open(output, "w") as f2:
        for line in f1:
            if line.startswith("#"):
                f2.write(line)
            else:
                elements = line.split()
                if (elements[0] == "DIOX_N" or elements[0] == "Dirigent") and elements[2].split("_")[2] == "lehmbachol":
                    f2.write(line)
                elif elements[0].startswith("Chal") and elements[2].split("_")[2] == "isogemichalcone":
                    f2.write(line)
                elif elements[0] in hmm_list and float(elements[4]) < 1e-50 and float(elements[7]) < 1e-50:
                    f2.write(line)

if __name__ == "__main__":
    hmmer_output = argv[1]
    hmmer_relevant_output = argv[2]

    select_relevant_hits(hmmer_output, hmmer_relevant_output)
