"""Gene-id harmonization — offline (mygene mocked)."""
import omics_oracle.genes as genes


def test_detect_id_type():
    assert genes.detect_id_type(["ENSG00000141510.5", "ENSG00000012048"]) == "ensembl"
    assert genes.detect_id_type(["TP53", "BRCA1", "EGFR"]) == "symbol"
    assert genes.detect_id_type(["7157", "672", "1956"]) == "entrez"
    assert genes.detect_id_type(["TP53", "ENSG00000012048", "672"]) == "mixed"


def test_strip_version():
    assert genes.strip_version("ENSG00000141510.5") == "ENSG00000141510"
    assert genes.strip_version("ENSG00000141510") == "ENSG00000141510"
    assert genes.strip_version("TP53") == "TP53"      # symbols untouched


def test_symbols_pass_through_without_network():
    mapping, ok = genes.map_to_symbols(["TP53", "brca1"], "symbol")
    assert ok and mapping["brca1"] == "BRCA1"


def test_ensembl_mapping_uses_mygene(monkeypatch, tmp_path):
    monkeypatch.setenv("OMICS_ORACLE_GENECACHE", str(tmp_path))
    # avoid reading a stale cache from a prior run
    monkeypatch.setattr(genes, "_load_cache", lambda species: {})
    def fake_query(ids, id_type, species):
        table = {"ENSG00000141510": "TP53", "ENSG00000012048": "BRCA1"}
        return {i: table.get(i, "") for i in ids}
    monkeypatch.setattr(genes, "_query_mygene", fake_query)
    mapping, ok = genes.map_to_symbols(["ENSG00000141510.5", "ENSG00000012048"], "ensembl")
    assert ok
    assert mapping["ENSG00000141510.5"] == "TP53"
    assert mapping["ENSG00000012048"] == "BRCA1"


def test_offline_falls_back_and_flags(monkeypatch, tmp_path):
    monkeypatch.setenv("OMICS_ORACLE_GENECACHE", str(tmp_path))
    monkeypatch.setattr(genes, "_load_cache", lambda species: {})
    monkeypatch.setattr(genes, "_query_mygene", lambda ids, t, s: None)  # network down
    mapping, ok = genes.map_to_symbols(["ENSG00000141510.5"], "ensembl")
    assert ok is False                                  # caller must warn
    assert mapping["ENSG00000141510.5"] == "ENSG00000141510"  # version-stripped fallback
