##############################################################################
#                Open Targets Rotation - Target Prioritisation               #
##############################################################################

# Packages
import os, pandas, json, pickle, numpy, argparse
pandas.options.mode.chained_assignment = None

##############################################################################
#                                I N P U T S                                 #
##############################################################################


def take_input():
    parser = argparse.ArgumentParser(prog="ProGoF - Prioritisation",
                                     usage="%(prog)s [inputs]",
                                     description="""
                                     **************************************
                                         Functional Consequence -> Drugs
                                             Drugs -> Prioritisation
                                     **************************************""")

    parser.add_argument("-version", dest="VERSION", required=True,
                        help="The interested Open Targets data version")

    parser.add_argument("-datap", dest="DATA_PATH", default=os.getcwd() + "/",
                        help="The path for input data.")

    parser.add_argument("-resultp", dest="RESULT_PATH", default=os.getcwd() + "/",
                        help="The path output data.")

    parsed_input = parser.parse_args()
    input_dict = vars(parsed_input)

    return input_dict


args = take_input()


##############################################################################
#                           F U N C T I O N S                                #
##############################################################################


def download_target_info(data_path, version):
    """
    Download whole target information in Open Targets
    :param data_path: The path in which the data will be downloaded
    :param version: Open Targets data version
    :return: True
    """

    download_url = "rsync -rpltvz --delete rsync.ebi.ac.uk::pub/databases/opentargets/" \
                   "platform/%s/output/etl/json/targets %s" % (version, data_path)

    os.system(download_url)

    return True


def read_target_info(data_path):
    """
    Read whole target information in Open Targets
    :param data_path: The path in which the data will be read
    :return: list of jsons
    """

    target_json_files = os.listdir(data_path + "targets/")
    jsons = []
    for json_file in target_json_files:
        json_str = open(data_path + "targets/" + json_file, 'r').read().replace("\n", " ")
        json_str = "[" + json_str + "]"
        json_str = json_str.replace("} {", "},{")
        json_obj = json.loads(json_str)
        for json_d in json_obj:
            jsons.append(json_d)

    return jsons


def add_target_info(data_path, version):
    """
    Construct the data frame from target information from Open Targets
    :param data_path: The path in which the data will be read
    :param version: Open Targets data version
    :return: Pandas Data Frame of the target information
    """

    if "targets" not in os.listdir(data_path):
        download_target_info(data_path=data_path, version=version)

    json_list = read_target_info(data_path=data_path)
    raw_json_df = pandas.json_normalize(json_list)

    raw_json_df.to_csv(data_path + "raw_target_df.csv")

    json_df = raw_json_df.copy()

    json_df["location"] = json_df.apply(lambda x: str(x["genomicLocation.chromosome"]) + ":" +
                                                  str(x["genomicLocation.start"]) + "-" +
                                                  str(x["genomicLocation.end"]), axis=1)

    json_df["antibody_buckets"] = json_df.apply(
        lambda x: ",".join([str(i) for i in x["tractability.antibody.buckets"]])
        if type(x["tractability.antibody.buckets"]) == list else "None", axis=1)

    json_df["small_molecule_buckets"] = json_df.apply(
        lambda x: ",".join([str(i) for i in x["tractability.smallmolecule.buckets"]])
        if type(x["tractability.smallmolecule.buckets"]) == list else "None", axis=1)

    json_df["subcellular_locations"] = json_df.apply(
        lambda x: ",".join(x["proteinAnnotations.subcellularLocations"])
        if type(x["proteinAnnotations.subcellularLocations"]) is not float else "None", axis=1)

    json_df["is_experimental_toxicity"] = json_df.apply(
        lambda x: True if type(x["safety.experimental_toxicity"]) is not float else False, axis=1)

    json_df["is_adverse_effects"] = json_df.apply(
        lambda x: True if type(x["safety.adverse_effects"]) is not float else False, axis=1)

    json_df = json_df[["id", "approvedSymbol", "bioType", "location", "genomicLocation.strand",
                       "subcellular_locations", "antibody_buckets", "small_molecule_buckets",
                       "tractability.antibody.categories.clinical_precedence",
                       "tractability.antibody.categories.predicted_tractable_high_confidence",
                       "tractability.antibody.categories.predicted_tractable_med_low_confidence",
                       "tractability.antibody.top_category",
                       "tractability.smallmolecule.categories.clinical_precedence",
                       "tractability.smallmolecule.categories.discovery_precedence",
                       "tractability.smallmolecule.categories.predicted_tractable",
                       "tractability.smallmolecule.top_category",
                       "is_experimental_toxicity", "is_adverse_effects"]]

    json_df.columns = ["ensembl_id", "gene_symbol", "bio_type", "location", "strand",
                       "subcellular_locations", "antibody_buckets", "small_molecule_buckets",
                       "antibody_clinical_precedence", "antibody_tractable_high_conf",
                       "antibody_tractable_med_low_conf", "antibody_category",
                       "small_molecule_clinical_precedence", "small_molecule_discovery_precedence",
                       "small_molecule_predicted_tractable", "small_molecule_category",
                       "is_experimental_toxicity", "is_adverse_effects"]

    json_df.to_csv(data_path + "curated_targets_df.csv")

    return json_df


def target_annotator(data_path, result_path, version):
    """
    Annotate the target set which were prioritised
    :param data_path: The path in which the data will be read
    :param result_path: The path in which the new annotated data will be written
    :param version: Open Targets data version
    :return: Pandas Data Frame of the target information
    """

    if "curated_targets_df.csv" not in os.listdir(data_path):
        target_df = add_target_info(data_path=data_path, version=version)
        target_df.to_csv(data_path + "curated_targets_df.csv")

    else:
        target_df = pandas.read_csv(data_path + "curated_targets_df.csv", index_col=0)

    target_df["min_tractability_bucket"] = "None"
    for i, row in target_df.iterrows():
        if row.antibody_buckets == "None" and row.small_molecule_buckets == "None":
            target_df.loc[i, "min_bucket"] = "None"
        else:
            if row.antibody_buckets == "None" and row.small_molecule_buckets != "None":
                if min([int(i) for i in row.small_molecule_buckets.split(",")]) in [1,2,3]:
                    target_df.loc[i, "min_tractability_bucket"] = "clinical"
                elif min([int(i) for i in row.small_molecule_buckets.split(",")]) in [4,5]:
                    target_df.loc[i, "min_tractability_bucket"] = "draggable"

            elif row.antibody_buckets != "None" and row.small_molecule_buckets == "None":
                if min([int(i) for i in row.antibody_buckets.split(",")]) in [1,2,3]:
                    target_df.loc[i, "min_tractability_bucket"] = "clinical"
                elif min([int(i) for i in row.antibody_buckets.split(",")]) in [4,5]:
                    target_df.loc[i, "min_tractability_bucket"] = "draggable"

            else:
                min_bucket = min([min([int(i) for i in row.antibody_buckets.split(",")]),
                                  min([int(i) for i in row.small_molecule_buckets.split(",")])])

                if min_bucket in [1,2,3]: target_df.loc[i, "min_tractability_bucket"] = "clinical"
                elif min_bucket in [4,5]: target_df.loc[i, "min_tractability_bucket"] = "draggable"

    target_df.to_csv(data_path + "curated_v2_targets_df.csv")

    validation_df = pandas.read_csv(result_path + "06_validation_merged_association_df.csv", index_col=0)
    validation_set = set(validation_df.gene_symbol)

    print("# targets in validation set : %s\n" % str(len(validation_set)))
    print("validation set : %s\n" % ",".join(list(validation_set)))

    merged_df = pandas.read_csv(result_path + "05_merged_association_df.csv", index_col=0)

    gene_prioritised_by_source = {}
    for i, row in merged_df.iterrows():
        gene = row.gene_symbol
        if gene not in gene_prioritised_by_source.keys():
            gene_prioritised_by_source[gene] = [row.datasource_id]
        else:
            if row.datasource_id not in gene_prioritised_by_source[gene]:
                gene_prioritised_by_source[gene].append(row.datasource_id)

    ellen_genes = pickle.load(open(args["DATA_PATH"] + "ellen_genes.txt", "rb"))
    ellen_gene_set = [gene for source, gene_list in ellen_genes.items() for gene in gene_list]

    df = pandas.DataFrame("None", index=set(gene_prioritised_by_source.keys()),
                          columns=["ensembl_id", "prioritised_sources", "is_Ellen",
                                   "in_validation", "tractability_by_min_tractability_bucket",
                                   "is_experimental_toxicity", "is_adverse_effects", "LoF_other"])

    for ind, _ in df.iterrows():
        m = merged_df[merged_df.gene_symbol == ind]
        df.loc[ind, "ensembl_id"] = list(set(m.ensembl_id))[0]
        x = target_df[target_df.ensembl_id == list(set(m.ensembl_id))[0]]
        df.loc[ind, "prioritised_sources"] = ",".join(gene_prioritised_by_source[ind])
        df.loc[ind, "is_Ellen"] = True if ind in ellen_gene_set else False
        df.loc[ind, "in_validation"] = True if ind in validation_set else False
        df.loc[ind, "tractability_by_min_tractability_bucket"] = list(x["min_tractability_bucket"])[0]
        df.loc[ind, "is_experimental_toxicity"] = list(x["is_experimental_toxicity"])[0]
        df.loc[ind, "is_adverse_effects"] = list(x["is_adverse_effects"])[0]
        df.loc[ind, "LoF_other"] = True if True in list(m["LoF_other"]) else False

    df.to_csv(result_path + "07_merged_target_df.csv")

    return df


def extract_gof_difference_lof(source, result_path):
    """
    Extract evidence having GoF/activating/increasing effect but
    not LoF/inhibiting/decreasing effect for the same disease
    :param source: The interested source type of the association
    :param result_path: The path in which the data will be written
    :return: The target list having interested effect but not opposite effect and
    the filtered association data frame
    """
    print("\n\n******************************************")
    print("%s GoF / LoF Filtration" % source.upper())
    print("******************************************\n")

    f = open(result_path + "stats_filtration_%s.txt" % source, "w")
    f.write("%s Statistics:\n" % source.upper())

    association_df = pandas.read_csv(result_path + "02_last_annotated_%s.csv" % source, index_col=0)

    print("Total target #: %d\n" % len(set(association_df.gene_symbol)))
    f.write("02 Association Data Frame:\n# Whole targets: %s\n" % str(len(set(association_df.gene_symbol))))

    gene_disease_pairs = set(association_df.apply(lambda x: x.gene_symbol + "-" + x.disease_id, axis=1))
    print("Total target-disease #: %d\n" % len(gene_disease_pairs))
    f.write("# Whole target-disease: %s\n" % str(len(gene_disease_pairs)))

    if source == "gene2phenotype": column_name = "SO_label"
    else: column_name = "functional_consq"

    whole_gof_annotations = [
        "Loss-of-function variants (as defined in pop up message) DO NOT cause this phenotype -please provide details in the comments",
        "Loss-of-function variants (as defined in pop up message) DO NOT cause this phenotype - please provide details in the comments",
        "copy_number_increase", "gain_of_function_variant", "protein gain of function",
        "gain_of_function_variant", "Increased function",
        "protein gain of function,variation affecting protein function"]

    whole_lof_annotations = ["effect on catalytic protein function,loss_of_function_variant",
                             "protein loss of function,unknown functional consequence",
                             "loss_of_function_variant", "decreased_translational_product_level",
                             "loss_of_function_variant,mutation affecting polypeptide function",
                             "loss_of_function_variant,sequence_variant_affecting_splicing",
                             "RNA degradation by nonsense-mediated decay,protein loss of function",
                             "Decreased function,Unknown function,variation affecting protein",
                             "Decreased function,variation affecting protein function",
                             "Decreased function", "effect on protein subcellular localization,protein truncation",
                             "protein loss of function", "complete gene deletion",
                             "polypeptide loss of function variant",
                             "RNA degradation by nonsense-mediated decay", "polypeptide_partial_loss_of_function",
                             "decreased_translational_product_level,effect on RNA splicing",
                             "loss_of_function_variant,mutation affecting polypeptide function,mutation affecting reading frame",
                             "amphigoric amino acid deletion", "loss_of_function_variant,protein loss of function",
                             "cryptic splice acceptor activation,protein loss of function", "protein truncation",
                             "protein truncation,variation affecting protein", "missing protein",
                             "decreased_translational_product_level,effect on protein degradation,loss_of_function_variant",
                             "decreased transcript level variant",
                             "RNA degradation by nonsense-mediated decay,effect on RNA splicing",
                             "reduced", "protein loss of function", "Decreased function", "loss_of_function_variant"]

    gof_association_df = association_df[association_df[column_name].isin(whole_gof_annotations)]
    print("# of Targets annotated as GoF/increase/activating: %d\n" % len(set(gof_association_df.gene_symbol)))
    f.write("# Whole targets with GoF/increase/activating: %s\n" % str(len(set(gof_association_df.gene_symbol))))

    gof_gene_diseases = {i: set(gof_association_df[gof_association_df.gene_symbol == i]["disease_id"])
                         for i in set(gof_association_df.gene_symbol)}

    pickle.dump(gof_gene_diseases, open(result_path + "gof_gene_diseases_%s.p" % source, "wb"))

    lof_association_df = association_df[association_df[column_name].isin(whole_lof_annotations)]
    print("# of Targets annotated as LoF/decrease/inhibiting: %d\n" % len(set(lof_association_df.gene_symbol)))
    f.write("# Whole targets with LoF/decrease/inhibiting: %s\n" % str(len(set(lof_association_df.gene_symbol))))

    lof_gene_diseases = {i: set(lof_association_df[lof_association_df.gene_symbol == i]["disease_id"])
                         for i in set(lof_association_df.gene_symbol)}

    pickle.dump(lof_gene_diseases, open(result_path + "lof_gene_diseases_%s.p" % source, "wb"))

    gof_not_lof_targets = []
    for target, disease_set in gof_gene_diseases.items():
        if target in lof_gene_diseases.keys():
            if len(disease_set.intersection(lof_gene_diseases[target])) == 0:
                gof_not_lof_targets.append(target)
        else:
            gof_not_lof_targets.append(target)

    print("# of Targets annotated as GoF/increase/activating but not \nLoF/decrease/inhibiting: %d\n"
          % len(set(gof_not_lof_targets)))
    f.write("# Whole targets with GoF but not LoF: %s\n" % str(len(set(gof_not_lof_targets))))
    f.close()
    target_filtered_df = association_df[(association_df.gene_symbol.isin(gof_not_lof_targets)) &
                                        (association_df[column_name] != "None") &
                                        (association_df[column_name].isin(whole_gof_annotations))]

    target_filtered_df.to_csv(result_path + "03_target_filtered_%s_df.csv" % source)

    return gof_not_lof_targets, target_filtered_df, gof_gene_diseases, lof_gene_diseases


def drug_mode_of_action(data_path):
    """
    Create a merged drug mode of action data frame
    :param data_path: The path in which the data will be read
    :return: Drug mode of action data frame including chembl_id of the drug,
    ensembl_id of the target, disease_id of the disease and
    mode of action of the drug on the corresponding disease and target
    """

    drug_df, moa_df = None, None

    if "core_drug_df.csv" not in os.listdir(data_path) or "drug_moa_df.csv" not in os.listdir(data_path):

        mof_json_files = os.listdir(data_path + "mechanismOfAction/")
        drug_json_files = os.listdir(data_path + "molecule/")

        mof_jsons, drug_jsons = [], []
        for mof_json_file in mof_json_files:
            mof_json_str = open(data_path + "mechanismOfAction/" + mof_json_file, 'r').read().replace("\n", " ")
            mof_json_str = "[" + mof_json_str + "]"
            mof_json_str = mof_json_str.replace("} {", "},{")
            mof_json_obj = json.loads(mof_json_str)
            for mof_json_d in mof_json_obj:
                mof_jsons.append(mof_json_d)

        for drug_json_file in drug_json_files:
            drug_json_str = open(data_path + "molecule/" + drug_json_file, 'r').read().replace("\n", " ")
            drug_json_str = "[" + drug_json_str + "]"
            drug_json_str = drug_json_str.replace("} {", "},{")
            drug_json_obj = json.loads(drug_json_str)
            for drug_json_d in drug_json_obj:
                drug_jsons.append(drug_json_d)

        raw_mof_json_df = pandas.json_normalize(mof_jsons)
        raw_drug_json_df = pandas.json_normalize(drug_jsons)

        # Core Drug information from Open Targets

        drug_json_df = raw_drug_json_df[["id", "drugType", "name", "hasBeenWithdrawn",
                                         "isApproved", "linkedTargets.rows", "linkedDiseases.rows"]]

        drug_json_df["linked_targets"] = drug_json_df.apply(
            lambda x: ",".join(x["linkedTargets.rows"]) if type(x["linkedTargets.rows"]) != float else (
                "None"), axis=1)

        drug_json_df["linked_diseases"] = drug_json_df.apply(
            lambda x: ",".join(x["linkedDiseases.rows"]) if type(x["linkedDiseases.rows"]) != float else (
                "None"), axis=1)

        drug_json_df = drug_json_df[["id", "drugType", "name", "hasBeenWithdrawn", "isApproved",
                                     "linked_targets", "linked_diseases"]]

        drug_json_df.columns = ["chembl_id", "drug_type", "drug_name", "has_been_withdrawn",
                                "is_approved", "linked_targets", "linked_diseases"]

        drug_json_df.to_csv(data_path + "core_drug_df.csv")

        # Drug mode of action information

        mof_json_df = raw_mof_json_df[["chemblIds", "targets", "targetType",
                                       "actionType", "mechanismOfAction", "references"]]

        moa_df = pandas.DataFrame(columns=["chembl_id", "linked_targets", "target_type",
                                           "action_type", "mechanism_of_action", "pmids"])

        for ind, row in mof_json_df.iterrows():
            targets = ",".join(row.targets) if type(row.targets) != float else "None"
            pmids = []
            for d in row.references:
                if d["source"] == "PubMed":
                    for pmid in d["ids"]:
                        pmids.append(pmid)

            if pmids: pmids_str = ",".join(pmids)
            else: pmids_str = "None"

            for chembl_id in row.chemblIds:
                moa_df = moa_df.append({"chembl_id": chembl_id,
                                        "linked_targets": targets,
                                        "target_type": row.targetType,
                                        "action_type": row.actionType,
                                        "mechanism_of_action": row.mechanismOfAction,
                                        "pmids": pmids_str}, ignore_index=True)

        moa_df.to_csv(data_path + "drug_moa_df.csv")

    else:
        drug_df = pandas.read_csv(data_path + "core_drug_df.csv", index_col=0)
        moa_df = pandas.read_csv(data_path + "drug_moa_df.csv", index_col=0)

    if "merged_moa_df.csv" not in os.listdir(data_path):
        merged_moa_df = pandas.DataFrame(columns = ["ensembl_id", "disease_id",
                                                    "chembl_id", "moa", "is_approved"])

        for ind, row in drug_df.iterrows():
            diseases = row["linked_diseases"].split(",")
            targets = row["linked_targets"].split(",")
            pairs = [(target, disease) for target in targets for disease in diseases]
            for pair in pairs:
                if pair[0] != "None" and pair[1] != "None":
                    merged_moa_df = merged_moa_df.append({"ensembl_id": pair[0], "disease_id": pair[1],
                                                          "chembl_id": row.chembl_id, "moa": "None",
                                                          "is_approved": row.is_approved},
                                                         ignore_index=True)

        for ind, row in moa_df.iterrows():

            if type(row["linked_targets"]) != float:
                targets = row["linked_targets"].split(",")
                for ensembl in targets:
                    indices = list(merged_moa_df[(merged_moa_df.ensembl_id == ensembl) &
                                                 (merged_moa_df.chembl_id == row.chembl_id)].index)

                    merged_moa_df.loc[indices, "moa"] = row.action_type

        merged_moa_df.to_csv(data_path + "merged_moa_df.csv")
    else:
        merged_moa_df = pandas.read_csv(data_path + "merged_moa_df.csv", index_col=0)

    return merged_moa_df


def annotate_tractability(data_path, result_path, source):
    """
    Annotate the lastly filtered association data frame with mode of action of the existing drugs
    :param data_path: The path in which the data will be read
    :param result_path: The path in which the data will be written
    :param source: The interested source type of the association
    :return:
    """

    # Read drug mode of action data frame
    merged_moa_df = drug_mode_of_action(data_path=data_path)

    # Read filtered annotated association data frame
    target_genes, target_filtered_df, gof, lof = extract_gof_difference_lof(source=source,
                                                                            result_path=result_path)

    # Create a new data frame since there can be several drugs and their mode of action
    # Write each of them separately
    moa_target_filtered_df = pandas.DataFrame(columns = list(target_filtered_df.columns) +
                                                        ["chembl_id", "moa", "is_approved"])
    moa_target_filtered_df["moa"] = "None"
    moa_target_filtered_df["chembl_id"] = "None"
    moa_target_filtered_df["is_approved"] = "None"

    index_counter = 0
    for ind, row in target_filtered_df.iterrows():
        # For each row in filtered annotated association data frame
        # Check if ensembl id and disease id pair is in drug mode of action data frame
        if (row.ensembl_id in set(merged_moa_df.ensembl_id)) and \
                (row.disease_id in set(merged_moa_df[merged_moa_df.ensembl_id == row.ensembl_id]["disease_id"])):
            for i, r in merged_moa_df[(merged_moa_df.ensembl_id == row.ensembl_id) &
                                      (merged_moa_df.disease_id == row.disease_id)].iterrows():
                # Each of the drug and mode of action pair for corresponding gene-disease pair
                # Write them on the new data frame
                for col in target_filtered_df.columns:
                    # Write first the existing columns
                    moa_target_filtered_df.loc[index_counter, col] = row[col]

                # Write, then the new columns
                moa_target_filtered_df.loc[index_counter, "chembl_id"] = r["chembl_id"]
                moa_target_filtered_df.loc[index_counter, "moa"] = r["moa"]
                moa_target_filtered_df.loc[index_counter, "is_approved"] = r["is_approved"]
                # Increment the index counter of the new data frame
                index_counter += 1

        else:
            # If there is no drug mode of action information, write None to that columns
            for col in target_filtered_df.columns:
                moa_target_filtered_df.loc[index_counter, col] = row[col]

            moa_target_filtered_df.loc[index_counter, "chembl_id"] = "None"
            moa_target_filtered_df.loc[index_counter, "moa"] = "None"
            moa_target_filtered_df.loc[index_counter, "is_approved"] = False

            index_counter += 1

    moa_target_filtered_df.to_csv(result_path + "04_drug_moa_%s.csv" % source)

    f = open(result_path + "stats_drug_%s.txt" % source, "w")
    f.write("%s Statistics:\n" % source.upper())

    # Targets having approved drugs having Chembl
    app_filtered_df = moa_target_filtered_df[(moa_target_filtered_df.is_approved) &
                                             (moa_target_filtered_df.chembl_id != "None")]
    print("# of targets having approved drugs having chembl id: %d\n"
          % len(set(app_filtered_df.gene_symbol)))
    f.write("# filtered targets having approved drugs having chembl id: %s\n" % str(len(set(app_filtered_df.gene_symbol))))

    # List the mode of actions that has negative/inhibitory effect on the target
    inhibitory_moas = ["NEGATIVE ALLOSTERIC MODULATOR", "BLOCKER", "INVERSE AGONIST", "DEGRADER",
                       "NEGATIVE MODULATOR", "INHIBITOR", "ANTAGONIST", "DISRUPTING AGENT", "ANTISENSE INHIBITOR"]

    # Filter the association data frame if they have approved inhibitory drug/molecule in Chembl
    inh_moa_app_filtered_df = app_filtered_df[app_filtered_df.moa.isin(inhibitory_moas)]

    print("# of targets having approved inhibitory chembl drugs/molecules: %d\n"
          % len(set(inh_moa_app_filtered_df.gene_symbol)))
    f.write("# filtered targets having approved inhibitory drugs having chembl id: %s\n" %
            str(len(set(inh_moa_app_filtered_df.gene_symbol))))
    f.close()

    return moa_target_filtered_df, target_filtered_df, gof, lof


def merging_prioritised_targets(data_path, result_path):

    merged_columns = ["ensembl_id", "gene_symbol", "disease_id", "disease_label", "score",
                      "datasource_id", "confidence", "variant_so_id", "is_approved",
                      "SO_label", "study_id", "variant_rsid", "functional_consq",
                      "last_clinvar_update", "panel_version", "panel_version_date",
                      "clinical_significances", "allelic_requirements", "chembl_id",
                      "moa", "tractability_type", "id"]

    if "source_df.p" not in os.listdir(result_path):
        source_info = {}
        for source in ["eva", "eva_somatic", "gene2phenotype", "genomics_england"]:

            moa_target_filtered_df, _, gof, lof = annotate_tractability(
                data_path=data_path, result_path=result_path, source=source)

            if source in ["eva", "eva_somatic"]:
                moa_target_filtered_df = moa_target_filtered_df.reset_index()
                if "index" in moa_target_filtered_df.columns:
                    moa_target_filtered_df = moa_target_filtered_df.rename(
                        columns={"index": "study_id"})

            source_df = moa_target_filtered_df[list(set(merged_columns).intersection(
                set(moa_target_filtered_df.columns)))]
            source_info[source] = {"df": source_df, "gof": gof, "lof": lof}

        pickle.dump(source_info, open(result_path + "source_df.p", "wb"))

    else:
        source_info = pickle.load(open(result_path + "source_df.p", "rb"))

    final_merged_df = pandas.DataFrame(
        index=range(sum([len(source_info[source]["df"].index) for source in source_info.keys()])),
        columns = merged_columns)

    index_counter = 0
    for source in source_info.keys():
        df = source_info[source]["df"]
        for ind, row in df.iterrows():
            final_merged_df.loc[index_counter, "ensembl_id"] = row.ensembl_id
            final_merged_df.loc[index_counter, "gene_symbol"] = row.gene_symbol
            final_merged_df.loc[index_counter, "disease_id"] = row.disease_id
            final_merged_df.loc[index_counter, "disease_label"] = row.disease_label
            final_merged_df.loc[index_counter, "score"] = row.score
            final_merged_df.loc[index_counter, "datasource_id"] = row.datasource_id
            final_merged_df.loc[index_counter, "confidence"] = row.confidence
            if source in ["eva", "eva_somatic", "gene2phenotype"]:
                final_merged_df.loc[index_counter, "variant_so_id"] = row.variant_so_id
                final_merged_df.loc[index_counter, "SO_label"] = row.SO_label
            else:
                final_merged_df.loc[index_counter, "variant_so_id"] = "None"
                final_merged_df.loc[index_counter, "SO_label"] = "None"

            final_merged_df.loc[index_counter, "study_id"] = row.study_id

            if source in ["eva", "eva_somatic"]:
                final_merged_df.loc[index_counter, "variant_rsid"] = row.variant_rsid
            else:
                final_merged_df.loc[index_counter, "variant_rsid"] = "None"

            if source in ["eva", "eva_somatic", "genomics_england"]:
                final_merged_df.loc[index_counter, "functional_consq"] = row.functional_consq
            else:
                final_merged_df.loc[index_counter, "functional_consq"] = "None"

            if source in ["eva", "eva_somatic"]:
                final_merged_df.loc[index_counter, "last_clinvar_update"] = row.last_clinvar_update
            else:
                final_merged_df.loc[index_counter, "last_clinvar_update"] = "None"

            if source in ["genomics_england"]:
                final_merged_df.loc[index_counter, "panel_version"] = row.panel_version \
                    if "panel_version" in df.columns else "None"
                final_merged_df.loc[index_counter, "panel_version_date"] = row.panel_version_date
            else:
                final_merged_df.loc[index_counter, "panel_version"] = "None"
                final_merged_df.loc[index_counter, "panel_version_date"] = "None"

            if source in ["eva", "eva_somatic"]:
                final_merged_df.loc[index_counter, "clinical_significances"] = row["clinvar_sig_(clinvar)"] \
                    if "clinvar_sig_(clinvar)" in df.columns else "None"
            else:
                final_merged_df.loc[index_counter, "clinical_significances"] = "None"

            if source in ["eva_somatic", "gene2phenotype"]:
                final_merged_df.loc[index_counter, "allelic_requirements"] = row.allelic_requirements \
                    if "allelic_requirements" in df.columns else "None"
            elif source in ["genomics_england"]:
                final_merged_df.loc[index_counter, "allelic_requirements"] = row.mode_of_inheritance \
                    if "mode_of_inheritance" in df.columns else "None"
            else:
                final_merged_df.loc[index_counter, "allelic_requirements"] = row.allele_origins \
                    if "allele_origins" in df.columns else "None"
            final_merged_df.loc[index_counter, "chembl_id"] = row.chembl_id
            final_merged_df.loc[index_counter, "moa"] = row.moa
            final_merged_df.loc[index_counter, "is_appoved"] = row.is_approved
            final_merged_df.loc[index_counter, "id"] = row.id
            index_counter += 1

    all_lof_pairs = []
    for source in source_info.keys():
        for gene, disease_list in source_info[source]["lof"].items():
            for disease in disease_list:
                if (gene, disease) not in all_lof_pairs:
                    all_lof_pairs.append((gene, disease))

    final_merged_df["LoF_other"] = "None"
    for ind, row in final_merged_df.iterrows():
        if (row.gene_symbol, row.disease_id) in all_lof_pairs:
            final_merged_df.loc[ind, "LoF_other"] = True
        else:
            final_merged_df.loc[ind, "LoF_other"] = False

    inhibitory_moas = ["NEGATIVE ALLOSTERIC MODULATOR", "BLOCKER", "INVERSE AGONIST", "DEGRADER",
                       "NEGATIVE MODULATOR", "INHIBITOR", "ANTAGONIST", "DISRUPTING AGENT", "ANTISENSE INHIBITOR"]

    final_merged_df["is_inhibitor"] = final_merged_df.apply(
        lambda x: True if x.moa in inhibitory_moas else False, axis=1)

    final_merged_df.to_csv(result_path + "05_merged_association_df.csv")

    validation_df = final_merged_df[final_merged_df.is_inhibitor]
    validation_df = validation_df[validation_df.LoF_other == False]

    validation_df.to_csv(result_path + "06_validation_merged_association_df.csv")

    return final_merged_df, validation_df


##############################################################################
#                      M A I N      E X E C U T I O N                        #
##############################################################################


def main(data_path, result_path, version):
    print("\n********************************\n")
    print("Prioritisation has been starting:\n")

    filtered, prioritised = merging_prioritised_targets(data_path=data_path,
                                                       result_path=result_path)

    target = target_annotator(data_path=data_path, result_path=result_path, version=version)

    print("Prioritisation has been finished.")
    print("********************************\n")

    return 1


main(data_path=args["DATA_PATH"], result_path=args["RESULT_PATH"], version = args["VERSION"])

