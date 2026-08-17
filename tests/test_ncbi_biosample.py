from http.client import IncompleteRead
from pathlib import Path

from genome_amr_predictor.data.ncbi_biosample import NCBIClient, parse_biosample_xml

FIXTURES = Path(__file__).parent / "fixtures"


def test_parse_structured_biosample_antibiogram():
    records = parse_biosample_xml((FIXTURES / "biosample_antibiogram.xml").read_bytes())

    assert len(records) == 6
    first = records[0]
    assert first.biosample_accession == "SAMN00000001"
    assert first.organism == "Escherichia coli"
    assert first.taxonomy_id == "562"
    assert first.bioproject_accession == "PRJNA000001"
    assert first.sra_accession == "SRS00000001"
    assert first.isolate == "EC-001"
    assert first.antibiotic == "ampicillin"
    assert first.phenotype_normalized == "resistant"
    assert first.binary_label == 1
    assert first.measurement == "32"
    assert first.testing_standard == "CLSI"


def test_parser_excludes_non_binary_outcomes_from_binary_label():
    records = parse_biosample_xml((FIXTURES / "biosample_antibiogram.xml").read_bytes())

    intermediate = next(record for record in records if record.antibiotic == "gentamicin")
    undefined = next(record for record in records if record.antibiotic == "meropenem")
    assert intermediate.phenotype_normalized == "intermediate"
    assert intermediate.binary_label is None
    assert undefined.phenotype_normalized == "excluded"
    assert undefined.binary_label is None


def test_large_eutils_request_uses_post(monkeypatch):
    captured = {}

    class Response:
        def __enter__(self):
            return self

        def __exit__(self, *_args):
            return False

        def read(self):
            return b"ok"

    def fake_urlopen(request, timeout):
        captured["request"] = request
        captured["timeout"] = timeout
        return Response()

    monkeypatch.setattr("urllib.request.urlopen", fake_urlopen)
    client = NCBIClient(timeout_seconds=5)

    payload = client._request("efetch.fcgi", {"id": ",".join(str(i) for i in range(500))})

    assert payload == b"ok"
    assert captured["request"].data is not None
    assert captured["request"].full_url.endswith("efetch.fcgi")
    assert captured["timeout"] == 5


def test_incomplete_response_is_retried(monkeypatch):
    attempts = 0

    class Response:
        def __enter__(self):
            return self

        def __exit__(self, *_args):
            return False

        def read(self):
            return b"complete"

    def flaky_urlopen(_request, timeout):
        nonlocal attempts
        assert timeout == 5
        attempts += 1
        if attempts == 1:
            raise IncompleteRead(b"partial")
        return Response()

    monkeypatch.setattr("urllib.request.urlopen", flaky_urlopen)
    monkeypatch.setattr("time.sleep", lambda _seconds: None)
    client = NCBIClient(timeout_seconds=5, retries=2)

    assert client._request("esearch.fcgi", {"term": "test"}) == b"complete"
    assert attempts == 2
