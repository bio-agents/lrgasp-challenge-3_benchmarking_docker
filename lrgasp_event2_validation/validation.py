#!/usr/bin/env python3
import os, io
import json
from argparse import ArgumentParser
from lrgasp import de_novo_rna_data
#from lrgasp import entry_metadata
#from lrgasp import experiment_metadata
from JSON_templates import JSON_templates

parser = ArgumentParser()
parser.add_argument("-i", "--input", help="input FASTA file for benchmarking", required=True)
parser.add_argument("-com", "--community_id", help="name of benchmarking community", required=True)
parser.add_argument("-c", "--challenges_ids", help="name of openEBench challenge", required=True)
parser.add_argument("-p", "--participant_id", help="name of the species used for prediction (mouse or manatee)",
                    required=True)
parser.add_argument("-r", "--public_ref_dir", help="directory with the list of genes used to validate the predictions",
                    required=True)
parser.add_argument("-o", "--output", help="output file where participant JSON file will be written",
                    required=True)
parser.add_argument("--coverage", help="identify possible splice-junctions", required=True)
args = parser.parse_args()


def main(args):
    # input parameters
    input = args.input
    public_ref_dir = args.public_ref_dir
    community_id = args.community_id
    challenges_ids = args.challenges_ids
    participant_id = args.participant_id
    splice_junctions = args.coverage + '/' + "*.tab"
    out_path = args.output

    # Assuring the output path does exist
    # if not os.path.exists(out_path):
    #    try:
    #        os.makedirs(out_path)
    #        with open(out_path, mode="a"):
    #            pass
    #    except OSError as exc:
    #       print("OS error: {0}".format(exc) + "\nCould not create output path: " + out_path)

    validate_input_data(input, public_ref_dir, community_id, challenges_ids, participant_id, out_path, splice_junctions)


def validate_input_data(input, public_ref_dir, community_id, challenges_ids, participant_id, out_path, splice_junctions):
    """
    Validates input data for the benchmarking community and the openEBench challenge.
    :param input: list of input files for benchmarking
    :param public_ref_dir: directory with the list of reference files
    :param coverage: identify possible splice-junctions
    :param community_id: name of benchmarking community
    :param challenges_ids: name of openEBench challenge
    :param participant_id: name of the species used for benchmarking
    :param out_path: output path where participant JSON file will be written
    :param experiment_json: Experiment JSON file that is requiered for uploading the submission. More info here: https://lrgasp.github.io/lrgasp-submissions/docs/metadata.html
    :param entry_json: Entry JSON file that is requiered for uploading the submission. More info here: https://lrgasp.github.io/lrgasp-submissions/docs/metadata.html
    """

    VALID_CHALLENGES_IDS = ['mouse_%_BUSCO_gene_found_vs_complete_CapTrap-ONT-LO',
           'mouse_%_BUSCO_gene_found_vs_complete_CapTrap-PacBio-LO', 'mouse_%_BUSCO_gene_found_vs_complete_R2C2-ONT-LO',
           'mouse_%_BUSCO_gene_found_vs_complete_cDNA-Illumina-SO', 'mouse_%_BUSCO_gene_found_vs_complete_cDNA-ONT-LO',
           'mouse_%_BUSCO_gene_found_vs_complete_cDNA-ONT-LS', 'mouse_%_BUSCO_gene_found_vs_complete_cDNA-PacBio-LO',
           'mouse_%_BUSCO_gene_found_vs_complete_cDNA-PacBio-LS', 'mouse_%_BUSCO_gene_found_vs_complete_dRNA-ONT-LO',
           'mouse_%_BUSCO_gene_found_vs_complete_dRNA-ONT-LS', 'mouse_%_CAGE_vs_%_QuantSeq_supp_trans_CapTrap-ONT-LO',
           'mouse_%_CAGE_vs_%_QuantSeq_supp_trans_CapTrap-PacBio-LO',
           'mouse_%_CAGE_vs_%_QuantSeq_supp_trans_R2C2-ONT-LO',
           'mouse_%_CAGE_vs_%_QuantSeq_supp_trans_cDNA-Illumina-SO',
           'mouse_%_CAGE_vs_%_QuantSeq_supp_trans_cDNA-ONT-LO', 'mouse_%_CAGE_vs_%_QuantSeq_supp_trans_cDNA-ONT-LS',
           'mouse_%_CAGE_vs_%_QuantSeq_supp_trans_cDNA-PacBio-LO',
           'mouse_%_CAGE_vs_%_QuantSeq_supp_trans_cDNA-PacBio-LS', 'mouse_%_CAGE_vs_%_QuantSeq_supp_trans_dRNA-ONT-LO',
           'mouse_%_CAGE_vs_%_QuantSeq_supp_trans_dRNA-ONT-LS',
           'mouse_%_canonical_SJ_vs_%_SJ_SR_coverage_CapTrap-ONT-LO',
           'mouse_%_canonical_SJ_vs_%_SJ_SR_coverage_CapTrap-PacBio-LO',
           'mouse_%_canonical_SJ_vs_%_SJ_SR_coverage_R2C2-ONT-LO',
           'mouse_%_canonical_SJ_vs_%_SJ_SR_coverage_cDNA-Illumina-SO',
           'mouse_%_canonical_SJ_vs_%_SJ_SR_coverage_cDNA-ONT-LO',
           'mouse_%_canonical_SJ_vs_%_SJ_SR_coverage_cDNA-ONT-LS',
           'mouse_%_canonical_SJ_vs_%_SJ_SR_coverage_cDNA-PacBio-LO',
           'mouse_%_canonical_SJ_vs_%_SJ_SR_coverage_cDNA-PacBio-LS',
           'mouse_%_canonical_SJ_vs_%_SJ_SR_coverage_dRNA-ONT-LO',
           'mouse_%_canonical_SJ_vs_%_SJ_SR_coverage_dRNA-ONT-LS',
           'mouse_%_trans_with_intra-priming_vs_rt_switching_CapTrap-ONT-LO',
           'mouse_%_trans_with_intra-priming_vs_rt_switching_CapTrap-PacBio-LO',
           'mouse_%_trans_with_intra-priming_vs_rt_switching_R2C2-ONT-LO',
           'mouse_%_trans_with_intra-priming_vs_rt_switching_cDNA-Illumina-SO',
           'mouse_%_trans_with_intra-priming_vs_rt_switching_cDNA-ONT-LO',
           'mouse_%_trans_with_intra-priming_vs_rt_switching_cDNA-ONT-LS',
           'mouse_%_trans_with_intra-priming_vs_rt_switching_cDNA-PacBio-LO',
           'mouse_%_trans_with_intra-priming_vs_rt_switching_cDNA-PacBio-LS',
           'mouse_%_trans_with_intra-priming_vs_rt_switching_dRNA-ONT-LO',
           'mouse_%_trans_with_intra-priming_vs_rt_switching_dRNA-ONT-LS', 'mouse_Antisense_CapTrap-ONT-LO',
           'mouse_Antisense_CapTrap-PacBio-LO', 'mouse_Antisense_R2C2-ONT-LO', 'mouse_Antisense_cDNA-Illumina-SO',
           'mouse_Antisense_cDNA-ONT-LO', 'mouse_Antisense_cDNA-ONT-LS', 'mouse_Antisense_cDNA-PacBio-LO',
           'mouse_Antisense_cDNA-PacBio-LS', 'mouse_Antisense_dRNA-ONT-LO', 'mouse_Antisense_dRNA-ONT-LS',
           'mouse_FSM_CapTrap-ONT-LO', 'mouse_FSM_CapTrap-PacBio-LO', 'mouse_FSM_R2C2-ONT-LO',
           'mouse_FSM_cDNA-Illumina-SO', 'mouse_FSM_cDNA-ONT-LO', 'mouse_FSM_cDNA-ONT-LS', 'mouse_FSM_cDNA-PacBio-LO',
           'mouse_FSM_cDNA-PacBio-LS', 'mouse_FSM_dRNA-ONT-LO', 'mouse_FSM_dRNA-ONT-LS', 'mouse_Fusion_CapTrap-ONT-LO',
           'mouse_Fusion_CapTrap-PacBio-LO', 'mouse_Fusion_R2C2-ONT-LO', 'mouse_Fusion_cDNA-Illumina-SO',
           'mouse_Fusion_cDNA-ONT-LO', 'mouse_Fusion_cDNA-ONT-LS', 'mouse_Fusion_cDNA-PacBio-LO',
           'mouse_Fusion_cDNA-PacBio-LS', 'mouse_Fusion_dRNA-ONT-LO', 'mouse_Fusion_dRNA-ONT-LS',
           'mouse_GenicGenomic_CapTrap-ONT-LO', 'mouse_GenicGenomic_CapTrap-PacBio-LO',
           'mouse_GenicGenomic_R2C2-ONT-LO', 'mouse_GenicGenomic_cDNA-Illumina-SO', 'mouse_GenicGenomic_cDNA-ONT-LO',
           'mouse_GenicGenomic_cDNA-ONT-LS', 'mouse_GenicGenomic_cDNA-PacBio-LO', 'mouse_GenicGenomic_cDNA-PacBio-LS',
           'mouse_GenicGenomic_dRNA-ONT-LO', 'mouse_GenicGenomic_dRNA-ONT-LS', 'mouse_GenicIntron_CapTrap-ONT-LO',
           'mouse_GenicIntron_CapTrap-PacBio-LO', 'mouse_GenicIntron_R2C2-ONT-LO', 'mouse_GenicIntron_cDNA-Illumina-SO',
           'mouse_GenicIntron_cDNA-ONT-LO', 'mouse_GenicIntron_cDNA-ONT-LS', 'mouse_GenicIntron_cDNA-PacBio-LO',
           'mouse_GenicIntron_cDNA-PacBio-LS', 'mouse_GenicIntron_dRNA-ONT-LO', 'mouse_GenicIntron_dRNA-ONT-LS',
           'mouse_ISM_CapTrap-ONT-LO', 'mouse_ISM_CapTrap-PacBio-LO', 'mouse_ISM_R2C2-ONT-LO',
           'mouse_ISM_cDNA-Illumina-SO', 'mouse_ISM_cDNA-ONT-LO', 'mouse_ISM_cDNA-ONT-LS', 'mouse_ISM_cDNA-PacBio-LO',
           'mouse_ISM_cDNA-PacBio-LS', 'mouse_ISM_dRNA-ONT-LO', 'mouse_ISM_dRNA-ONT-LS',
           'mouse_Intergenic_CapTrap-ONT-LO', 'mouse_Intergenic_CapTrap-PacBio-LO', 'mouse_Intergenic_R2C2-ONT-LO',
           'mouse_Intergenic_cDNA-Illumina-SO', 'mouse_Intergenic_cDNA-ONT-LO', 'mouse_Intergenic_cDNA-ONT-LS',
           'mouse_Intergenic_cDNA-PacBio-LO', 'mouse_Intergenic_cDNA-PacBio-LS', 'mouse_Intergenic_dRNA-ONT-LO',
           'mouse_Intergenic_dRNA-ONT-LS', 'mouse_NIC_CapTrap-ONT-LO', 'mouse_NIC_CapTrap-PacBio-LO',
           'mouse_NIC_R2C2-ONT-LO', 'mouse_NIC_cDNA-Illumina-SO', 'mouse_NIC_cDNA-ONT-LO', 'mouse_NIC_cDNA-ONT-LS',
           'mouse_NIC_cDNA-PacBio-LO', 'mouse_NIC_cDNA-PacBio-LS', 'mouse_NIC_dRNA-ONT-LO', 'mouse_NIC_dRNA-ONT-LS',
           'mouse_NNC_CapTrap-ONT-LO', 'mouse_NNC_CapTrap-PacBio-LO', 'mouse_NNC_R2C2-ONT-LO',
           'mouse_NNC_cDNA-Illumina-SO', 'mouse_NNC_cDNA-ONT-LO', 'mouse_NNC_cDNA-ONT-LS', 'mouse_NNC_cDNA-PacBio-LO',
           'mouse_NNC_cDNA-PacBio-LS', 'mouse_NNC_dRNA-ONT-LO', 'mouse_NNC_dRNA-ONT-LS',
           'mouse_num_trans_vs_with_coding_potential_CapTrap-ONT-LO',
           'mouse_num_trans_vs_with_coding_potential_CapTrap-PacBio-LO',
           'mouse_num_trans_vs_with_coding_potential_R2C2-ONT-LO',
           'mouse_num_trans_vs_with_coding_potential_cDNA-Illumina-SO',
           'mouse_num_trans_vs_with_coding_potential_cDNA-ONT-LO',
           'mouse_num_trans_vs_with_coding_potential_cDNA-ONT-LS',
           'mouse_num_trans_vs_with_coding_potential_cDNA-PacBio-LO',
           'mouse_num_trans_vs_with_coding_potential_cDNA-PacBio-LS',
           'mouse_num_trans_vs_with_coding_potential_dRNA-ONT-LO',
           'mouse_num_trans_vs_with_coding_potential_dRNA-ONT-LS', 'mouse_sirvs_CapTrap-ONT-LO',
           'mouse_sirvs_CapTrap-PacBio-LO', 'mouse_sirvs_R2C2-ONT-LO', 'mouse_sirvs_cDNA-Illumina-SO',
           'mouse_sirvs_cDNA-ONT-LO', 'mouse_sirvs_cDNA-ONT-LS', 'mouse_sirvs_cDNA-PacBio-LO',
           'mouse_sirvs_cDNA-PacBio-LS', 'mouse_sirvs_dRNA-ONT-LO', 'mouse_sirvs_dRNA-ONT-LS',
           'mouse_trans_CapTrap-ONT-LO', 'mouse_trans_CapTrap-PacBio-LO', 'mouse_trans_R2C2-ONT-LO',
           'mouse_trans_cDNA-Illumina-SO', 'mouse_trans_cDNA-ONT-LO', 'mouse_trans_cDNA-ONT-LS',
           'mouse_trans_cDNA-PacBio-LO', 'mouse_trans_cDNA-PacBio-LS', 'mouse_trans_dRNA-ONT-LO',
           'mouse_trans_dRNA-ONT-LS']

    # validate challenges ids
    challenges_ids_split = challenges_ids.split(" ")
    for challenge_id in challenges_ids_split:
        if challenge_id not in VALID_CHALLENGES_IDS:
            print("Challenge ID not valid: " + challenge_id)
            exit(1)

    # validate if input file is FASTA format
    validated = False

    try:
        de_novo_rna_data.load(input)
    except Exception as e:
        print("Error: " + str(e) + "\nInput file is not in FASTA format.")
        exit(1)

    # check if public reference data exists
    public_ref_annotation = public_ref_dir + '/' + "lrgasp_gencode_vM27_sirvs.gtf"
    public_ref_genome = public_ref_dir + '/' + "lrgasp_grcm39_sirvs.fasta"

    if not os.path.exists(public_ref_annotation):
        print("Error: " + public_ref_annotation + " does not exist.")
        exit(1)

    if not os.path.exists(public_ref_genome):
        print("Error: " + public_ref_genome + " does not exist.")
        exit(1)

    # check if splice-junctions file exists
    if not os.path.exists(splice_junctions):
        print("Warning: " + splice_junctions + " does not exist.")

    # check if experiment JSON file exists
    #if not os.path.exists(experiment_json):
    #    print("Warning: " + "experiment_json file:" + experiment_json + " does not exist, use a fake one.")
    #else:
    #    try:
    #        experiment_metadata.load(experiment_json)
    #    except Exception as e:
    #        print("Error: " + str(e) + "\nExperiment JSON file is not in LRGASP JSON format.More info here: "
    #                                   "https://lrgasp.github.io/lrgasp-submissions/docs/metadata.html")

    # check if entry JSON file exists
    #if not os.path.exists(entry_json):
    #    print("Warning: " + "entry_json file:" + entry_json + " does not exist, use a fake one.")
    #    ID = "fake_ID"
    #else:
    #    try:
    #        entry_data = entry_metadata.load(entry_json)
    #        ID = entry_data['entry_id']
    #    except Exception as e:
    #        print("Error: " + str(e) + "\nEntry JSON file is not in LRGASP JSON format.More info here: "
    #                                   "https://lrgasp.github.io/lrgasp-submissions/docs/metadata.html")
    #        ID = "fake_ID"

    validated = True
    ID = "fake_ID"
    # write participant JSON file
    if validated:
        print("Validation Finished! Input and metadata is valid. Validation result is written to " + out_path)
        output_json = JSON_templates.write_participant_dataset(ID, community_id, challenges_ids, participant_id, validated)        # write participant JSON file
        with open(out_path, 'w') as outfile:
            json.dump(output_json, outfile, sort_keys=True, indent=4, separators=(',', ': '))


if __name__ == "__main__":
    main(args)
