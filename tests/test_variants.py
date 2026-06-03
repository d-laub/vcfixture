from vcfixture import variants as v
from vcfixture.allele import Seq


def test_constructors():
    assert v.snp("A", "T") == ("A", "T")
    assert v.mnp("AC", "GT") == ("AC", "GT")
    assert v.insertion("A", "CGT") == ("A", "ACGT")
    assert v.deletion("A", "CGT") == ("ACGT", "A")
    assert v.delins("AC", "GTA") == ("AC", "GTA")
    assert v.spanning_deletion("A") == ("A", "*")


def test_classify_single_alt():
    assert v.classify("A", "T") == "SNP"
    assert v.classify("AC", "GT") == "MNP"
    assert v.classify("A", "ACGT") == "INS"
    assert v.classify("ACGT", "A") == "DEL"
    assert v.classify("AC", "GTA") == "DELINS"
    assert v.classify("A", "*") == "SPANNING_DEL"


def test_record_class_multiallelic():
    assert v.record_class("G", (Seq("A"), Seq("C"))) == "MULTIALLELIC"
    assert v.record_class("G", (Seq("A"),)) == "SNP"
