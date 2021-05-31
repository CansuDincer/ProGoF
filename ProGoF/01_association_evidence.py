##############################################################################
#                  Open Targets Rotation - Data download                     #
##############################################################################

# Packages
import os, pandas, json, requests, pickle, argparse, numpy
import xml.etree.ElementTree as ET

pandas.options.mode.chained_assignment = None


##############################################################################
#                               I N P U T S                                  #
##############################################################################


def take_input():
    parser = argparse.ArgumentParser(prog="ProGof - Retrieval & Annotation",
                                     usage="%(prog)s [inputs]",
                                     description="""
                                     **************************************
                                     Open Targets | Gene->Disease->Evidence
                                        Evidence->Functional Consequence
                                     **************************************""")

    parser.add_argument("-source", dest="SOURCE", required=True,
                        choices=["genomics_england", "gene2phenotype",
                                 "eva_somatic", "eva", "ot_genetics_portal"],
                        help="The source of the target-disease association.")

    parser.add_argument("-version", dest="VERSION", required=True,
                        help="The interested Open Targets data version")

    parser.add_argument("-score", dest="SCORE", required=True,
                        help="Threshold for the genetic association score."
                             "The given number and greater of them will be taken.")

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


def download_associations(source, data_path, version):
    """
    Download direct associations from Open Targets by source type
    :param source: The interested source type of the association
    :param data_path: The path in which the data will be downloaded
    :param version: Open Targets data version
    :return: True
    """

    download_url = "rsync -rpltvz --delete rsync.ebi.ac.uk::pub/databases/opentargets/platform/" \
                   "%s/output/etl/json/evidence/sourceId=%s %s" % (version, source, data_path)

    os.system(download_url)

    return True


def read_associations(source, data_path):
    """
    Read direct associations from Open Targets by source type
    :param source: The interested source type of the association
    :param data_path: The path in which the data will be read
    :return: list of jsons
    """

    association_json_files = os.listdir(data_path + "sourceID=%s/" % source)
    jsons = []
    for json_file in association_json_files:
        json_str = open(data_path + "sourceID=" + source + "/" + json_file, 'r').read().replace("\n", " ")
        json_str = "[" + json_str + "]"
        json_str = json_str.replace("} {", "},{")
        json_obj = json.loads(json_str)
        for json_d in json_obj:
            jsons.append(json_d)

    return jsons


def association_dataframes(source, data_path, version):
    """
    Construct the data frame from direct associations from Open Targets
    :param source: The interested source type of the association
    :param data_path: The path in which the data will be read
    :param version: Open Targets data version
    :return: Pandas Data Frame of the associations for interested source
    """

    if "sourceId=%s" % source not in os.listdir(data_path):
        download_associations(source=source, data_path=data_path, version=version)

    json_list = read_associations(source=source, data_path=data_path)
    raw_json_df = pandas.json_normalize(json_list)

    raw_json_df.to_csv(data_path + "raw_%s_df.csv" % source)

    if source == "gene2phenotype":
        json_df = raw_json_df.copy()
        json_df["allelic_requirements"] = json_df.apply(lambda x: ",".join(
            x["allelicRequirements"]) if type(x["allelicRequirements"]) is list else x["allelicRequirements"], axis=1)

        json_df["pmids"] = json_df.apply(lambda x: ",".join(
            x["literature"]) if type(x["literature"]) is list else x["literature"], axis=1)

        json_df = json_df[["targetId", "targetSymbol", "variantFunctionalConsequenceId",
                           "diseaseId", "diseaseLabel", "score", "datasourceId", "confidence",
                           "allelic_requirements", "studyId", "id", "datatypeId", "pmids"]]

        json_df.columns = ["ensembl_id", "gene_symbol", "variant_so_id",
                           "disease_id", "disease_label", "score", "datasource_id", "confidence",
                           "allelic_requirements", "study_id", "id", "datatype_id", "pmids"]

        return json_df

    elif source == "eva":
        json_df = raw_json_df.copy()
        json_df["allele_origins"] = json_df.apply(
            lambda x: ",".join(x["alleleOrigins"]) if type(x["alleleOrigins"]) is list
            else x["alleleOrigins"], axis=1)

        json_df["clinical_significances"] = json_df.apply(
            lambda x: ",".join(x["clinicalSignificances"]) if type(x["clinicalSignificances"]) is list
            else x["clinicalSignificances"], axis=1)

        json_df["pmids"] = json_df.apply(
            lambda x: ",".join(x["literature"]) if type(x["literature"]) is list
            else x["literature"], axis=1)

        json_df = json_df[["targetId", "targetSymbol", "variantFunctionalConsequenceId",
                           "diseaseId", "diseaseLabel", "score", "datasourceId", "confidence",
                           "variantRsId", "variantId", "clinical_significances",
                           "allele_origins", "studyId", "id", "datatypeId", "pmids"]]

        json_df.columns = ["ensembl_id", "gene_symbol", "variant_so_id",
                           "disease_id", "disease_label", "score", "datasource_id", "confidence",
                           "variant_rsid", "variant_id", "clinical_significances",
                           "allele_origins", "study_id", "id", "datatype_id", "pmids"]

        return json_df

    elif source == "eva_somatic":
        json_df = raw_json_df.copy()
        json_df["allele_origins"] = json_df.apply(
            lambda x: ",".join(x["alleleOrigins"]) if type(x["alleleOrigins"]) is list else x["alleleOrigins"], axis=1)
        json_df["clinical_significances"] = json_df.apply(
            lambda x: ",".join(x["clinicalSignificances"]) if type(x["clinicalSignificances"]) is list else x[
                "clinicalSignificances"], axis=1)
        json_df["pmids"] = json_df.apply(
            lambda x: ",".join(x["literature"]) if type(x["literature"]) is list else x["literature"], axis=1)

        json_df["allelic_requirements"] = json_df.apply(lambda x: ",".join(
            x["allelicRequirements"]) if type(x["allelicRequirements"]) is list else x["allelicRequirements"], axis=1)

        json_df = json_df[["targetId", "targetSymbol", "variantFunctionalConsequenceId", "diseaseId",
                           "diseaseLabel", "score", "datasourceId", "confidence", "variantRsId",
                           "variantId", "clinical_significances", "allele_origins",
                           "allelic_requirements", "studyId", "id", "datatypeId", "pmids"]]

        json_df.columns = ["ensembl_id", "gene_symbol", "variant_so_id", "disease_id",
                           "disease_label", "score", "datasource_id", "confidence", "variant_rsid",
                           "variant_id", "clinical_significances", "allele_origins",
                           "allelic_requirements", "study_id", "id", "datatype_id", "pmids"]

        return json_df

    elif source == "genomics_england":
        json_df = raw_json_df.copy()
        json_df["pmids"] = json_df.apply(
            lambda x: ",".join(x["literature"]) if type(x["literature"]) is list else x["literature"], axis=1)

        json_df["allelic_requirements"] = json_df.apply(lambda x: ",".join(
            x["allelicRequirements"]) if type(x["allelicRequirements"]) is list else x["allelicRequirements"], axis=1)

        json_df = json_df[["targetId", "targetSymbol", "diseaseId", "diseaseLabel", "score",
                           "datasourceId", "confidence", "allelic_requirements", "studyId",
                           "id", "datatypeId", "pmids"]]

        json_df.columns = ["ensembl_id", "gene_symbol", "disease_id", "disease_label", "score",
                           "datasource_id", "confidence", "allelic_requirements", "study_id",
                           "id", "datatype_id", "pmids"]

        return json_df

    elif source == "ot_genetics_portal":
        json_df = raw_json_df.copy()
        json_df["pmids"] = json_df.apply(
            lambda x: ",".join(x["literature"]) if type(x["literature"]) is list else x["literature"], axis=1)

        json_df = json_df[["targetId", "targetSymbol", "variantFunctionalConsequenceId",
                           "diseaseId", "diseaseLabel", "score", "studySampleSize", "beta",
                           "oddsRatio", "pValueExponent", "pValueMantissa", "datasourceId",
                           "variantRsId", "variantId", "studyId", "id",
                           "datatypeId", "publicationFirstAuthor", "publicationYear", "pmids"]]

        json_df.columns = ["ensembl_id", "gene_symbol", "variant_so_id", "disease_id",
                           "disease_label", "score", "study_sample_size", "beta", "odds_ratio",
                           "exponent_p", "Mantissa_p", "datasource_id", "variant_rsid",
                           "variant_id", "study_id", "id", "datatype_id", "pubs_first_author",
                           "pubs_year", "pmids"]

        return json_df

    else:
        return False


def so_ontology_retriever(data_path):
    """
    Read Sequence Ontology (SO) information.
    :param data_path: The path in which the data will be read
    :return: Python dictionary for ontology information
    """

    ontologies = {}
    so_ont_str = open(data_path + "ontology-so-2021-04-16.json", "r").read().replace("\n", " ")
    so_ont_str = "[" + so_ont_str + "]"
    so_ont_str = so_ont_str.replace("} {", "},{")
    so_ont_obj = json.loads(so_ont_str)
    for json_d in so_ont_obj:
        so_id = json_d["id"].replace(":", "_")
        if so_id not in ontologies.keys():
            ontologies[so_id] = [json_d]
        else:
            t = ontologies[so_id]
            if json_d not in t:
                t.append(json_d)

            ontologies[so_id] = t

    return ontologies


def so_ontology_analyser(df, ontology_dict):
    """
    Annotate the evidence of the association data frame with Sequence Ontology
    :param df: Association Data Frame
    :param ontology_dict: The python dictionary for ontology information - so_ontology_retriever()
    :return: Pandas Data Frame of the SO information
    """

    df["SO_label"] = "None"
    for so_val in list(set(df["variant_so_id"])):
        if so_val in ontology_dict.keys():
            if "label" in ontology_dict[so_val][0].keys():
                indices = list(df[df["variant_so_id"] == so_val].index)
                label = ontology_dict[so_val][0]["label"]
                df.loc[indices, "SO_label"] = label

    return df


def clinvar_rcv_retriever(rcv):
    """
    Collection of variant information through ClinVar API
    :param rcv: Variant Disease Record id
    :return: XML root (element tree format)
    """

    clinvar_rcv_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" \
                      "db=clinvar&rettype=clinvarset&id="

    request = clinvar_rcv_url + rcv

    response = requests.get(request)

    if response.status_code != 200:
        print("No response from ensembl!")
        return 0

    else:
        root = ET.fromstring(response.content)
        print(rcv)
        return root


def clinvar_rcv_analyser(root):
    """
    Python style storing of variant information from ClinVar XML
    :param root: XML root (element tree format)
    :return: Python dictionary having variant information
    """
    clinvar_info = {}

    for child in root:

        # Latest update time
        clinvar_info["last_update"] = child.find("ReferenceClinVarAssertion").attrib["DateLastUpdated"]

        # RCV accession to check with our RCVs
        clinvar_info["accession"] = child.find("ReferenceClinVarAssertion/ClinVarAccession").attrib["Acc"]

        # Clinical significance of the variant
        clinvar_info["clinvar_sig"] = child.find("ReferenceClinVarAssertion/ClinicalSignificance/"
                                                 "Description").text

        # The source of the significance and any annotations
        clinvar_info["source_type"] = child.find("ReferenceClinVarAssertion/ObservedIn/Method/MethodType").text

        # Functional and molecular consequences of the variants if any
        cons_node = child.find("ReferenceClinVarAssertion/MeasureSet/Measure").findall("AttributeSet")
        fcons_list, mcons_list = [], []
        for node in cons_node:
            # Functional Consequences
            cons_node2 = node.find("Attribute[@Type='FunctionalConsequence']")
            if cons_node2 is not None:
                fcons_list.append(cons_node2.text)
            else:
                continue

            # Molecular Consequences
            cons_node3 = node.find("Attribute[@Type='MolecularConsequence']")
            if cons_node3 is not None:
                mcons_list.apend(cons_node3.text)
            else:
                continue

        if not fcons_list:
            clinvar_info["functional_consq"] = "None"
        else:
            t = ",".join(fcons_list)
            if t[-1] == ",": t = t[:-1]
            clinvar_info["functional_consq"] = t

        if not mcons_list:
            clinvar_info["molecular_consq"] = "None"
        else:
            t = ",".join(mcons_list)
            if t[-1] == ",": t = t[:-1]
            clinvar_info["molecular_consq"] = t

        # Disease mechanism of the variants
        trait_node = child.find("ReferenceClinVarAssertion/TraitSet/Trait").findall("AttributeSet")
        trait_list = []

        for node in trait_node:
            # Disease mechanism
            trait_node2 = node.find("Attribute[@Type='disease mechanism']")
            if trait_node2 is not None:
                trait_list.append(trait_node2.text)
            else:
                continue

        if not trait_list:
            clinvar_info["disease_mech"] = "None"
        else:
            t = ",".join(trait_list)
            if t[-1] == ",": t = t[:-1]
            clinvar_info["disease_mech"] = t

    return clinvar_info


def xref_rcv_rsid(rsid_list, data_path, result_path, result_xref_dict):
    """
    Obtaining RCV ids of ClinVar of corresponding rs ids of dbSNP
    :param rsid_list: List of rs ids from association data frame
    :param data_path: Path of the input data (Xref file)
    :param result_path: Path of the output data - python dictionary
    :param result_xref_dict: Name of the resulting file
    :return: Python dictionary-pickle format - key: rs id - value: list of RCV id(s)
    """

    xref_df = pandas.read_csv(data_path + "variant_summary.txt", sep="\t")
    xref_df["rsid"] = xref_df.apply(lambda x: "rs" + str(x["RS# (dbSNP)"]), axis=1)
    list_of_variants_dict = {}
    for rs in rsid_list:
        if rs in list(xref_df["rsid"]):
            if rs not in list_of_variants_dict.keys():
                list_of_variants_dict[rs] = {}
                for gene, rcv_list in xref_df[xref_df.rsid == rs].groupby(["GeneSymbol", "RCVaccession"]).groups:
                    rcv_total_list = []
                    for rcv in rcv_list.split("|"):
                        if rcv not in rcv_total_list:
                            rcv_total_list.append(rcv)
                    if gene not in list_of_variants_dict[rs].keys():
                        list_of_variants_dict[rs][gene] = rcv_total_list
                    else:
                        k = list_of_variants_dict[rs][gene]
                        for rcv in rcv_total_list:
                            if rcv not in k:
                                k.append(rcv)
                        list_of_variants_dict[rs][gene] = k

            else:
                t = list_of_variants_dict[rs]
                for gene, rcv_list in xref_df[xref_df.rsid == rs].groupby(["GeneSymbol", "RCVaccession"]).groups:
                    if gene in t.keys():
                        for rcv in rcv_list.split("|"):
                            if rcv not in t[gene]:
                                t[gene].append(rcv)
                    else:
                        t[gene] = {}
                        rcv_total_list = []
                        for rcv in rcv_list.split("|"):
                            if rcv not in rcv_total_list:
                                rcv_total_list.append(rcv)
                        t[gene] = rcv_total_list

                list_of_variants_dict[rs] = t

    pickle.dump(list_of_variants_dict, open(result_path + result_xref_dict + ".p", "wb"))

    return list_of_variants_dict


def clinvar_annotator(df, variant_type, xref_dict, data_path, result_path, source):
    """
    Annotate the evidence of the association data frame with ClinVar
    :param df: Association Data Frame annotated with SO annotation
    :param variant_type: The type of the variant identifier - rsid or RCV
    :param xref_dict: Python dictionary-pickle format - key: rs id - value: list of RCV id(s) - xref_rcv_rsid()
    :param data_path: Path of the rcv data
    :param result_path: The path in which the new SO annotated data will be written
    :param source: The interested source type of the association
    :return: Pandas Data Frame of the SO and ClinVar annotation
    """

    new_df, rcv_column = None, None

    if variant_type == "RCV":
        rcv_column = "study_id"
        new_df = df.copy()

    elif variant_type == "rsid":
        variant_column = "variant_rsid"
        rcv_df = pandas.DataFrame()
        for rsid in list(set(df[variant_column])):
            df_row = df[df[variant_column] == rsid]
            df_row["RCV"] = "None"
            if rsid in xref_dict.keys():
                for gene, rcv_list in xref_dict[rsid].items():
                    if gene in list(set(df_row["gene_symbol"])):
                        gene_df_row = df_row[df_row.gene_symbol == gene]
                        gene_df_row = gene_df_row.append([gene_df_row] * (len(rcv_list) - 1), ignore_index=True)
                        for k in range(0, len(rcv_list)):
                            gene_df_row.loc[k, "RCV"] = rcv_list[k]
                        rcv_df = pandas.concat([rcv_df, gene_df_row])

        genes_yes_rcv = list(set(rcv_df["gene_symbol"]))
        genes_no_rcv_df = df[~df["gene_symbol"].isin(genes_yes_rcv)]
        genes_no_rcv_df["RCV"] = "None"
        new_df = pandas.concat([rcv_df, genes_no_rcv_df])
        rcv_column = "RCV"

    new_df = new_df.set_index(rcv_column)
    new_df["last_clinvar_update"] = "None"
    new_df["clinvar_sig_(clinvar)"] = "None"
    new_df["source_type"] = "None"
    new_df["functional_consq"] = "None"
    new_df["molecular_consq"] = "None"
    new_df["disease_mech"] = "None"
    new_df["VCF_Allele"] = "None"
    clinvar_dict = {}

    if "rcvs" not in os.listdir(data_path):
        os.system("mkdir rcvs %s." % result_path)

    k_total = len(set(list(new_df.index)))
    k = 0
    # Downloading the RCV files since this process os taking so long
    # For not to start from the beginning
    # Additionally, showing the percentage to follow
    for rcv in set(list(new_df.index)):
        k += 1
        if rcv != "None":
            if "%s_%s_dict.p" % (rcv, source) not in os.listdir(result_path + "rcvs/"):
                root = clinvar_rcv_retriever(rcv=rcv)
                if root:
                    d = clinvar_rcv_analyser(root)
                    if d != {}:
                        clinvar_dict[rcv] = d
                        pickle.dump(d, open(result_path + "rcvs/%s_%s_dict.p" % (rcv, source), "wb"))
            else:
                d = pickle.load(open(result_path + "rcvs/%s_%s_dict.p" % (rcv, source), "rb"))
                clinvar_dict[rcv] = d

        print(k * 100.0 / k_total)

    pickle.dump(clinvar_dict, open(result_path + "ClinVar_RCV_%s_dict.p" % source, "wb"))

    clinvar_annotated_df = pandas.DataFrame()
    for rcv, var_dict in clinvar_dict.items():
        raw_df = new_df[new_df.index == rcv]
        raw_df["last_clinvar_update"] = var_dict["last_update"]
        raw_df["clinvar_sig_(clinvar)"] = var_dict["clinvar_sig"]
        raw_df["source_type"] = var_dict["source_type"]
        raw_df["functional_consq"] = var_dict["functional_consq"]
        raw_df["molecular_consq"] = var_dict["molecular_consq"]
        raw_df["disease_mech"] = var_dict["disease_mech"]
        raw_df["VCF_Allele"] = var_dict["VCF_Allele"] \
            if "VCF_Allele" in var_dict.keys() else "None"
        clinvar_annotated_df = pandas.concat([clinvar_annotated_df, raw_df])

    genes_yes_clinvar = list(set(clinvar_annotated_df["gene_symbol"]))
    genes_no_clinvar_df = new_df[~new_df["gene_symbol"].isin(genes_yes_clinvar)]
    new_clinvar_df = pandas.concat([clinvar_annotated_df, genes_no_clinvar_df])

    return new_clinvar_df


def panelapp_retriever(gene_symbol):
    """
    The API retrieval of the interested gene's information from PanelAPP.
    :param gene_symbol: Hugo symbol of the interested gene
    :return: json response of the interested gene (only elements under results key)
    """

    panelapp_gene_server = "https://panelapp.genomicsengland.co.uk/api/v1/genes/"

    response = requests.get(panelapp_gene_server + gene_symbol + "/",
                            headers={"Content-Type": "application/json"})

    if response.status_code != 200:
        print("No response from PanelAPP!")
        return 0

    else:
        return response.json()["results"]


def panelapp_analyser(gene_symbol, panel):
    """
    The analysis of the PanelAPP retrieval-mode of pathogenicity
    :param gene_symbol: Hugo symbol of the interested gene
    :param panel: The interested panel id
    :return: Python dictionary including all gene-panel related information
    """

    panelapp_info = dict()
    result_list = panelapp_retriever(gene_symbol=gene_symbol)
    for result_dict in result_list:
        panel_id = result_dict["panel"]["id"]
        if int(panel_id) == int(panel):
            panelapp_info["confidence_level"] = result_dict["confidence_level"]
            panelapp_info["functional_consq"] = result_dict["mode_of_pathogenicity"]
            panelapp_info["mode_of_inheritance"] = result_dict["mode_of_inheritance"]
            panelapp_info["panel_version"] = result_dict["panel"]["version"]
            panelapp_info["panel_version_date"] = result_dict["panel"]["version_created"]
            panelapp_info["panel_pmids"] = ",".join(result_dict["publications"])

    return panelapp_info


def panelapp_annotator(df):
    """
    Annotation of the Genomics England dataframe with PanelApp API
    :param df: Genomics England association data frame
    :return: New panel app annotated Genomics England association data frame
    """

    new_df = df.copy()

    new_df["confidence_level"] = "None"
    new_df["functional_consq"] = "None"
    new_df["mode_of_inheritance"] = "None"
    new_df["functional_consq"] = "None"
    new_df["panel_version"] = "None"
    new_df["panel_version_date"] = "None"
    new_df["panel_pmids"] = "None"

    for ind, row in df.iterrows():
        d = panelapp_analyser(gene_symbol=row.gene_symbol, panel=row.study_id)
        if d != {}:
            new_df.at[ind, "confidence_level"] = d["confidence_level"]
            new_df.at[ind, "functional_consq"] = d["functional_consq"]
            new_df.at[ind, "mode_of_inheritance"] = d["mode_of_inheritance"]
            new_df.at[ind, "functional_consq"] = d["functional_consq"]
            new_df.at[ind, "panel_version"] = d["panel_version"]
            new_df.at[ind, "panel_version_date"] = d["panel_version_date"]
            new_df.at[ind, "panel_pmids"] = d["panel_pmids"]

    return new_df


##############################################################################
#                      M A I N      E X E C U T I O N                        #
##############################################################################


def main(source, data_path, result_path, version, score):
    print("\n************************")
    print("%s" % source.upper())
    print("************************\n")
    print("Analysis has been starting:\n")

    f = open(result_path + "stats_%s.txt" % source, "w")
    f.write("%s Statistics:\n" % source.upper())

    # Association data frames

    if "curated_%s_df.csv" % source not in os.listdir(result_path):

        df = association_dataframes(source=source, data_path=data_path, version=version)

    else:
        df = pandas.read_csv(result_path + "01_curated_%s_df.csv" % source, index_col=0)

    whole_targets = len(set(df.gene_symbol))
    f.write("01 Curated Data Frame:\n# Whole targets: %s\n" % str(whole_targets))

    whole_diseases = len(set(df.disease_id))
    f.write("# Whole diseases: %s\n" % str(whole_diseases))

    whole_associations = len(set([(row.gene_symbol, row.disease_id) for i, row in df.iterrows()]))
    f.write("# Whole associations: %s\n" % str(whole_associations))

    print("Association data frame was created.\n")

    print("Filtering out the association below the score as %s." % score)

    df = df[df.score >= float(score)]

    score_filtered_targets = len(set(df.gene_symbol))
    f.write("01 Curated Data Frame after score filtration:\n# Whole targets: %s\n" % str(score_filtered_targets))

    score_filtered_diseases = len(set(df.disease_id))
    f.write("# Whole diseases: %s\n" % str(score_filtered_diseases))

    score_filtered_associations = len(set([(row.gene_symbol, row.disease_id) for i, row in df.iterrows()]))
    f.write("# Whole associations: %s\n" % str(score_filtered_associations))

    print("\nAnnotation is starting:")

    if source in ["ot_genetics_portal", "eva", "eva_somatic", "gene2phenotype"]:

        print("\nSO annotation:")

        ontology_dict = so_ontology_retriever(data_path=data_path)

        so_annotated_df = so_ontology_analyser(df=df, ontology_dict=ontology_dict)

        so_annotated_targets = len(set(so_annotated_df[so_annotated_df.SO_label != "None"].gene_symbol))
        f.write("\nSequence Ontology annotation:\n")
        f.write("# SO annotated targets: %s\n" % str(so_annotated_targets))

        so_annotated_diseases = len(set(so_annotated_df[so_annotated_df.SO_label != "None"].disease_id))
        f.write("# SO annotated diseases: %s\n" % str(so_annotated_diseases))

        so_annotated_associations = len(set(
            [(row.gene_symbol, row.disease_id)
             for i, row in so_annotated_df[so_annotated_df.SO_label != "None"].iterrows()]))
        f.write("# SO annotated associations: %s\n" % str(so_annotated_associations))

        print("     SO annotation was added.\n")

        if source == "ot_genetics_portal":

            print("\nClinVar annotation:")

            # OT Genetics Portal do not have RCV id but dbSNP RS id
            # There is a function to collect corresponding RCV id for ClinVar analysis
            print("OT Genetics RS ids will be converted to RCV ids for analysis.")

            xref_dict = xref_rcv_rsid(rsid_list=set(so_annotated_df.variant_rsid),
                                      data_path=data_path, result_path=result_path,
                                      result_xref_dict="ot_genetics_rsid_rcv")

            last_df = clinvar_annotator(df=so_annotated_df, variant_type="rsid", data_path=data_path,
                                        xref_dict=xref_dict, result_path=result_path, source=source)

            last_df.to_csv(result_path + "02_last_annotated_%s.csv" % source)

            clinvar_targets = len(set(last_df[last_df.functional_consq != "None"].gene_symbol))
            f.write("\nClinVar annotation:\n")
            f.write("# ClinVar annotated targets: %s\n" % str(clinvar_targets))

            clinvar_diseases = len(set(last_df[last_df.functional_consq != "None"].disease_id))
            f.write("# ClinVar annotated diseases: %s\n" % str(clinvar_diseases))

            clinvar_associations = len(set(
                [(row.gene_symbol, row.disease_id)
                 for i, row in last_df[last_df.functional_consq != "None"].iterrows()]))
            f.write("# ClinVar annotated associations: %s\n" % str(clinvar_associations))

            print("     ClinVar annotation was added.\n")

        elif source in ["eva", "eva_somatic"]:

            print("\nClinVar annotation:")

            # We directly have RCV ids of the corresponding variants --> directly into analysis
            last_df = clinvar_annotator(df=so_annotated_df, variant_type="RCV", data_path=data_path,
                                        xref_dict=None, result_path=result_path, source=source)

            last_df.to_csv(result_path + "02_last_annotated_%s.csv" % source)

            clinvar_targets = len(set(last_df[last_df.functional_consq != "None"].gene_symbol))
            f.write("\nClinVar annotation:\n")
            f.write("# ClinVar annotated targets: %s\n" % str(clinvar_targets))

            clinvar_diseases = len(set(last_df[last_df.functional_consq != "None"].disease_id))
            f.write("# ClinVar annotated diseases: %s\n" % str(clinvar_diseases))

            clinvar_associations = len(set(
                [(row.gene_symbol, row.disease_id)
                 for i, row in last_df[last_df.functional_consq != "None"].iterrows()]))
            f.write("# ClinVar annotated associations: %s\n" % str(clinvar_associations))

            print("     ClinVar annotation was added.\n")

        elif source == "gene2phenotype":

            last_df = so_annotated_df.copy()

            last_df.to_csv(result_path + "02_last_annotated_%s.csv" % source)

    if source == "genomics_england":

        print("\nPanelApp annotation:")

        last_df = panelapp_annotator(df=df)

        last_df.to_csv(result_path + "02_last_annotated_%s.csv" % source)

        panelapp_targets = len(set(last_df[last_df.functional_consq != "None"].gene_symbol))
        f.write("\nPanelApp annotation:\n")
        f.write("# PanelApp annotated targets: %s\n" % str(panelapp_targets))

        panelapp_diseases = len(set(last_df[last_df.functional_consq != "None"].disease_id))
        f.write("# PanelApp annotated diseases: %s\n" % str(panelapp_diseases))

        panelapp_associations = len(set(
            [(row.gene_symbol, row.disease_id)
             for i, row in last_df[last_df.functional_consq != "None"].iterrows()]))
        f.write("# PanelApp annotated associations: %s\n" % str(panelapp_associations))

        print("     PanelApp annotation was added.\n")

    print("Analysis has been finished for %s!" % source.upper())
    print("************************\n")

    return True


main(source=args["SOURCE"], data_path=args["DATA_PATH"], result_path=args["RESULT_PATH"],
     version=args["VERSION"], score=args["SCORE"])


