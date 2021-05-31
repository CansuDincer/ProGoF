##############################################################################
#             Open Targets Rotation - Benchmarking and analysis              #
##############################################################################

# Packages
import os, pandas, json, pickle, numpy
from matplotlib import pyplot as plt
from matplotlib_venn import venn2, venn2_unweighted, venn2_circles
from matplotlib_venn import venn3, venn3_unweighted, venn3_circles
import matplotlib.patches as patches
import venn

args = {"RESULT_PATH": "/Users/cd7/Documents/Rotations/Trynka_Group/opentargets/venv/data/deneme2/output/",
        "DATA_PATH": "/Users/cd7/Documents/Rotations/Trynka_Group/opentargets/venv/data/deneme2/input/"}

sources = ["genomics_england", "gene2phenotype","eva_somatic", "eva"]

##############################################################################
#                      A N N O T A T I O N   S T A T S                       #
##############################################################################


##############################################################################
#                      F I N A L    D A T A F R A M E S                      #
##############################################################################

final_merged_df = pandas.read_csv(args["RESULT_PATH"] + "05_merged_association_df.csv")

final_merged_target_df = pandas.read_csv(args["RESULT_PATH"] + "07_merged_target_df.csv",index_col = 0)

# All raw targets
all_raw_targets = []
all_raw_diseases = []
all_raw_associations = []

for source in sources:
    raw_df = pandas.read_csv(args["DATA_PATH"] + "raw_%s_df.csv" % source, index_col=0)
    genes = set(raw_df.targetSymbol)
    for gene in genes:
        if gene not in all_raw_targets: all_raw_targets.append(gene)

    diseases = set(raw_df.diseaseId)
    for disease in diseases:
        if disease not in all_raw_diseases: all_raw_diseases.append(disease)

    associations = set(raw_df.groupby(["targetSymbol", "diseaseId"]).groups.keys())
    for association in associations:
        if association not in all_raw_associations: all_raw_associations.append(association)


# All targets
all_targets = []
all_diseases = []
all_associations = []

for source in sources:
    curated_df = pandas.read_csv(args["RESULT_PATH"] + "02_last_annotated_%s.csv" % source,
                                 index_col=0)
    genes = set(curated_df.gene_symbol)
    for gene in genes:
        if gene not in all_targets: all_targets.append(gene)

    diseases = set(curated_df.disease_id)
    for disease in diseases:
        if disease not in all_diseases: all_diseases.append(disease)

    associations = set(curated_df.groupby(["gene_symbol", "disease_id"]).groups.keys())
    for association in associations:
        if association not in all_associations: all_associations.append(association)


for source in sources:
    raw_df = pandas.read_csv(args["DATA_PATH"] + "raw_%s_df.csv" % source, index_col=0)
    curated_df = pandas.read_csv(args["RESULT_PATH"] + "02_last_annotated_%s.csv" % source,
                                 index_col=0)
    curated_df = curated_df.reset_index()

    filtered_df = pandas.read_csv(args["RESULT_PATH"] +
                                  "03_target_filtered_%s_df.csv" % source, index_col=0)

    if source != "gene2phenotype":
        evidence_id = "study_id"
        consequence = "functional_consq"
    else:
        evidence_id = "variant_so_id"
        consequence = "SO_label"

    curated_df[consequence] = curated_df[consequence].fillna("None")

    print(source)
    print("Raw targets: %d" % len(set(raw_df.targetSymbol)))
    print("Raw diseases: %d" % len(set(raw_df.diseaseId)))
    print("Raw associations: %d" % len(raw_df.groupby(["targetSymbol", "diseaseId"]).size()))
    print("0.4 filtered targets: %d" % len(set(curated_df.gene_symbol)))
    print("0.4 filtered diseases: %d" % len(set(curated_df.disease_id)))
    print("0.4 filtered associations: %d" % len(curated_df.groupby(["gene_symbol", "disease_id"]).size()))
    annotated_df = curated_df[curated_df[consequence] != "None"]
    print("0.4 filtered-annotated targets: %d" % len(set(annotated_df.gene_symbol)))
    print("0.4 filtered-annotated diseases: %d" % len(set(annotated_df.disease_id)))
    print("0.4 filtered-annotated associations: %d" % len(
        annotated_df.groupby(["gene_symbol", "disease_id"]).size()))

    print("0.4 filtered-annotated association evidences: %d" % len(
        annotated_df.groupby(["gene_symbol", "disease_id", evidence_id]).size()))
    print("Unique evidence identifier: %d" % len(set(annotated_df[evidence_id])))

    print("GoF filtered targets: %d" % len(set(filtered_df.gene_symbol)))
    print("GoF filtered diseases: %d" % len(set(filtered_df.disease_id)))
    print("GoF filtered associations: %d" % len(
        filtered_df.groupby(["gene_symbol", "disease_id"]).size()))


# Filtered targets / GoF \ LoF
filtered_targets = set(final_merged_df.gene_symbol)

# Validation targets from prioritised targets
# They have inhibitory drugs with a Chembl ID
validation_targets = set(final_merged_target_df[final_merged_target_df.in_validation].index)

# All prioritised targets /
# draggable and no experimental toxicity and adverse effects and no LoF from another source
all_prioritised = set(final_merged_target_df[(~final_merged_target_df.is_adverse_effects) &
                                             (~final_merged_target_df.is_experimental_toxicity) &
                                             (~final_merged_target_df.LoF_other) &
                                             (final_merged_target_df.tractability_by_min_tractability_bucket != "None")].index)

non_draggable = filtered_targets.difference(all_prioritised.union(validation_targets))

filtered_targets_disvby_sources = final_merged_target_df.groupby(["prioritised_sources"]).size()
prioritised_targets_disvby_sources = final_merged_target_df[final_merged_target_df.index.isin(all_prioritised)]\
    .groupby(["prioritised_sources"]).size()


##############################################################################
#        B E N C H M A R K I N G   w / E L L E N ' S    R E S U L T S        #
##############################################################################

ellen_genes = pickle.load(open(args["DATA_PATH"] + "ellen_genes.txt", "rb"))
all_ellen_targets = set([target for source, target_list in ellen_genes.items() for target in target_list])

# Ellen vs all targets
ellen_and_all_targets = all_ellen_targets.intersection(set(all_targets))
ellen_not_all_targets = all_ellen_targets.difference(set(all_targets))

# Ellen vs filtered targets
ellen_and_filtered = all_ellen_targets.intersection(filtered_targets)
ellen_not_filtered = all_ellen_targets.difference(filtered_targets)
filtered_not_ellen = filtered_targets.difference(all_ellen_targets)

# Ellen vs prioritised targets
ellen_and_prioritised = all_ellen_targets.intersection(all_prioritised)
ellen_not_prioritised = all_ellen_targets.difference(all_prioritised)
prioritised_not_ellen = all_prioritised.difference(all_ellen_targets)


def draw_venn3(labels, names, title, savefile, fs):

    v3_uw = venn3_unweighted(labels, set_labels=names, set_colors=('w', 'w', 'w'))
    v3_c = venn3_circles(subsets=(1, 1, 1, 1, 1, 1, 1), linewidth=2, color="black")

    for text in v3_uw.set_labels:
        text.set_fontsize(fs)
        text.set_fontname("Calibri")
    for text in v3_uw.subset_labels:
        text.set_fontsize(fs)
        text.set_fontname("Calibri")
    plt.title(title, fontname="Calibri", fontsize=20)
    plt.tight_layout()
    plt.savefig(args["RESULT_PATH"] + "figures/%s.pdf" % savefile, dpi=300)
    plt.savefig(args["RESULT_PATH"] + "figures/%s.jpg" % savefile, dpi=300)
    plt.close()
    return 1

# Venn for Ellen comparison

# All - Ellen - Filtered
labels = [all_ellen_targets, set(all_targets), filtered_targets]

names = ["Ellen's", "ProGoF initial targets", "Filtered"]

draw_venn3(labels=labels, names=names, title="Comparison of number of Targets",
           savefile="Venn_Ellen_All_Filtered", fs = 16)

# Ellen - Filtered - Prioritised
labels = [all_ellen_targets, filtered_targets, all_prioritised]

names = ["Ellen's", "Filtered", "Prioritised"]

draw_venn3(labels=labels, names=names, title="Comparison of number of Targets",
           savefile="Venn_Ellen_Filtered_Prioritised", fs = 16)

# Filtered : Prioritised - Druggable - None

labels = [validation_targets, all_prioritised.difference(validation_targets),
          filtered_targets.difference(all_prioritised)]

names = ["Validation", "Drug Development Promise", "No Drug"]

draw_venn3(labels=labels, names=names, title="GoF Filtered Targets",
           savefile="Venn_inside_filtered", fs = 16)


# Stats

curated_ge_df = pandas.read_csv(args["RESULT_PATH"] + "02_last_annotated_%s.csv" % "genomics_england",
                                index_col=0)
curated_g2p_df = pandas.read_csv(args["RESULT_PATH"] + "02_last_annotated_%s.csv" % "gene2phenotype",
                                 index_col=0)

filtered_ge_df = pandas.read_csv(args["RESULT_PATH"] + "03_target_filtered_%s_df.csv" % "genomics_england",
                                 index_col=0)
filtered_g2p_df = pandas.read_csv(args["RESULT_PATH"] + "03_target_filtered_%s_df.csv" % "gene2phenotype",
                                  index_col=0)


# Intersection and differences

g2p_ellen_and_progof = set(ellen_genes["Gene2Phenotype"]).intersection(set(curated_g2p_df.gene_symbol))
g2p_ellen_not_progof = set(ellen_genes["Gene2Phenotype"]).difference(set(curated_g2p_df.gene_symbol))
g2p_progof_not_ellen = set(curated_g2p_df.gene_symbol).difference(set(ellen_genes["Gene2Phenotype"]))
filter_g2p_ellen_and_progof = set(ellen_genes["Gene2Phenotype"]).intersection(set(filtered_g2p_df.gene_symbol))
filter_g2p_ellen_not_progof = set(ellen_genes["Gene2Phenotype"]).difference(set(filtered_g2p_df.gene_symbol))
filter_g2p_ellen_not_progof.difference(g2p_ellen_not_progof)
filter_g2p_progof_not_ellen = set(filtered_g2p_df.gene_symbol).difference(set(ellen_genes["Gene2Phenotype"]))

ge_ellen_and_progof = set(ellen_genes["PanelApp"]).intersection(set(curated_ge_df.gene_symbol))
ge_ellen_not_progof = set(ellen_genes["PanelApp"]).difference(set(curated_ge_df.gene_symbol))
ge_progof_not_ellen = set(curated_ge_df.gene_symbol).difference(set(ellen_genes["PanelApp"]))
filter_ge_ellen_and_progof = set(ellen_genes["PanelApp"]).intersection(set(filtered_ge_df.gene_symbol))
filter_ge_ellen_not_progof = set(ellen_genes["PanelApp"]).difference(set(filtered_ge_df.gene_symbol))
filter_ge_ellen_not_progof.difference(ge_ellen_not_progof)
filter_ge_progof_not_ellen = set(filtered_ge_df.gene_symbol).difference(set(ellen_genes["PanelApp"]))


all_ellen_targets.difference(set(curated_g2p_df.gene_symbol).union(set(curated_ge_df.gene_symbol)))
all_ellen_targets.difference(all_ellen_targets.difference(set(curated_g2p_df.gene_symbol).union(set(curated_ge_df.gene_symbol))).union(set(filtered_ge_df.gene_symbol).union(set(filtered_g2p_df.gene_symbol))))
