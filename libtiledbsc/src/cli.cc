// This file is currently a sandbox for C++ API experiments

#include <tiledbsc/tiledbsc>

using namespace tiledbsc;

void walk_soco(std::string_view uri) {
    auto soco = SOMACollection::open(uri);
    auto somas = soco->list_somas();

    LOG_INFO(fmt::format("walking soco URI = {}", uri));

    for (auto& [name, uri] : somas) {
        LOG_INFO(fmt::format("  soma {} = {}", name, uri));

        auto soma = SOMA::open(uri);
        auto arrays = soma->list_arrays();
        for (auto& [name, uri] : arrays) {
            LOG_INFO(fmt::format("    array {} = {}", name, uri));
        }
    }
}

void slice_soma(std::string_view soma_uri) {
    Config conf;
    conf["config.logging_level"] = "5";

    auto soma = SOMA::open(soma_uri, conf);
    auto array = soma->open_array("obs");
    auto mq = ManagedQuery(array);

    auto q = std::make_unique<Query>(array->schema().context(), *array);
    auto est_bytes = q->est_result_size_var("obs_id");
    LOG_DEBUG(fmt::format("est_num_cells = {}", est_bytes[0] / 8));
    LOG_DEBUG(fmt::format("est_bytes = {}", est_bytes[1]));

    mq.select_columns({"obs_id", "percent_mito"});
    mq.select_points<std::string>("obs_id", {"AAACATACAACCAC-1"});
    mq.select_ranges<std::string>(
        "obs_id", {{"TTTCGAACTCTCAT-1", "TTTGCATGCCTCAC-1"}});

    size_t total_cells = 0;
    while (!mq.is_complete()) {
        auto num_cells = mq.submit();
        LOG_DEBUG(fmt::format("num_cells = {}", num_cells));
        total_cells += num_cells;
        continue;

        auto mito = mq.data<float>("percent_mito");
        for (size_t i = 0; i < num_cells; i++) {
            auto obs = mq.string_view("obs_id", i);
            LOG_INFO(
                fmt::format("obs_id = {} percent_mito = {}", obs, mito[i]));
        }
    }
    LOG_DEBUG(fmt::format("total_cells = {}", total_cells));
}

void soco_query(std::string_view soco_uri) {
    Config config;
    // config["config.logging_level"] = "5";
    // config["soma.init_buffer_bytes"] = "4294967296";
    // config["soma.init_buffer_bytes"] = "268435456";

    auto soco = SOMACollection::open(soco_uri, config);
    auto sqs = soco->query();
    auto ctx = soco->context();

    std::vector<std::string> obs_attrs = {
        "assay",
        "assay_ontology_term_id",
        "cell_type",
        "cell_type_ontology_term_id",
        "development_stage",
        "development_stage_ontology_term_id",
        "disease",
        "disease_ontology_term_id",
        "ethnicity",
        "ethnicity_ontology_term_id",
        "is_primary_data",
        "organism",
        "organism_ontology_term_id",
        "sex",
        "sex_ontology_term_id",
        "tissue",
        "tissue_ontology_term_id"};

    std::vector<std::string> var_attrs = {
        "feature_biotype",
        "feature_is_filtered",
        "feature_name",
        "feature_reference"};

    std::string obs_attr = "cell_type";
    std::string obs_val = "pericyte cell";
    auto obs_qc = QueryCondition::create(*ctx, obs_attr, obs_val, TILEDB_EQ);
    sqs->set_obs_condition(obs_qc);
    sqs->select_obs_attrs(obs_attrs);

    std::string var_attr = "feature_name";
    std::string var_val = "DPM1";
    auto var_qc = QueryCondition::create(*ctx, var_attr, var_val, TILEDB_EQ);
    // sqs->set_var_condition(var_qc);
    sqs->select_var_attrs(var_attrs);

    while (auto results = sqs->next_results()) {
    }

    LOG_DEBUG("Done");
}

int main(int argc, char** argv) {
    if (argc != 3) {
        printf("Usage: %s soco_uri soma_uri\n", argv[0]);
        return 1;
    }

    LOG_CONFIG("debug");

    try {
        // walk_soco(argv[1]);
        // slice_soma(argv[2]);
        // debug(argv[2]);
        soco_query(argv[1]);
    } catch (const std::exception& e) {
        LOG_FATAL(e.what());
    }

    return 0;
};
