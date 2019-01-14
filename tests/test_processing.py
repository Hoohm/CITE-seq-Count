import pytest
from cite_seq_count import processing


@pytest.fixture
def data():
    import json
    from collections import defaultdict
    from collections import OrderedDict
    from collections import Counter
    
    from itertools import islice
    # Test file paths
    pytest.correct_R1_path = 'tests/test_data/fastq/correct_R1.fastq.gz'
    pytest.correct_R2_path = 'tests/test_data/fastq/correct_R2.fastq.gz'
    
    pytest.chunk_size = 800
    pytest.tags = OrderedDict({
        'CGTCGTAGCTGATCGTAGCTGAC': 'test3-CGTCGTAGCTGATCGTAGCTGAC',
        'CGTACGTAGCCTAGC': 'test2-CGTACGTAGCCTAGC',
        'CGTAGCTCG': 'test1-CGTAGCTCG'
        })
    pytest.barcode_slice = slice(0, 16)
    pytest.umi_slice = slice(16, 26)
    pytest.first_line = 1
    pytest.correct_whitelist = set(['ACTGTTTTATTGGCCT','TTCATAAGGTAGGGAT'])
    pytest.legacy = False
    pytest.debug = False
    pytest.start_trim = 0
    pytest.maximum_distance = 5
    pytest.corrected_results = {
        'ACTGTTTTATTGGCCT':
        defaultdict(None, {
            'test3-CGTCGTAGCTGATCGTAGCTGAC': Counter({b'CATTAGTGGT': 2, b'ACACTCNTGT': 2, b'AAGCCTTTAN': 1, b'CCCACGGGTG': 1, b'AAACAAAAGA': 1, b'ACTTCTTAGC': 1, b'AGAGNCCTGA': 1, b'GACNATTGCA': 1, b'TTCAATNGGC': 1, b'CNCCTCGCTC': 1, b'AACGTCGTTT': 1, b'CNGNTATCTT': 1, b'GCTATTCTGG': 1, b'CCAGATTAGT': 1, b'GNTGTCCTGA': 1, b'TCTAAACNAN': 1, b'NAANNTTATA': 1, b'CCCCGGCAAT': 1, b'GATAGTNCGA': 1, b'TATGGTNTTA': 1, b'GTGGTGATAT': 1, b'GATTATCTTC': 1, b'TCGCTTGTAT': 1, b'CAGAACTTCT': 1, b'ACAGGCCACA': 1, b'TTACTTCCTC': 1, b'ATCTCTCAAA': 1, b'CACGTTANAG': 1, b'GCCCANATTA': 1, b'CGGGAGTGCG': 1}),
            'test1-CGTAGCTCG': Counter({b'CNCCTCGCTC': 2, b'CTGTGGGACA': 2, b'AAAATCGGNG': 2, b'AAGCCTTTAN': 2, b'TCTTNATCCT': 1, b'GTGGTGATAT': 1, b'TTGAATACTG': 1, b'ANTACATGTA': 1, b'GANCAGCCTC': 1, b'ATAGGGTTGC': 1, b'ATAACCCCTC': 1, b'CAGAACTTCT': 1, b'AGTATAATAC': 1, b'TATGCTTTGG': 1, b'TGAAATCTTG': 1, b'GGGTTAGCGC': 1, b'ACTTCTTAGC': 1, b'CTTTTNCCNG': 1, b'AATTCAAATC': 1, b'AGAGNCCTGA': 1, b'CATTAGTGGT': 1, b'GNTGTCCTGA': 1, b'GCTACAAGGG': 1, b'TAGTGCTGAG': 1, b'TCGGTCACAA': 1, b'AAACGCGCCG': 1, b'GAAANCATTA': 1, b'CCTACAAATA': 1, b'TTGCAGTGAT': 1, b'CGGGAGTGCG': 1, b'AAGTTTCTCT': 1, b'CGAGGGGCTG': 1}),
            'test2-CGTACGTAGCCTAGC':
            Counter({b'TCGTATGGCN': 2, b'CGGGAGTGCG': 2, b'CCAGGTGGGA': 2, b'TTACTTCCTC': 2, b'GANCAGCCTC': 1, b'CATTAGTGGT': 1, b'AATTCAAATC': 1, b'CACGTTANAG': 1, b'CNGNTATCTT': 1, b'NAANNTTATA': 1, b'TCANCCNCCG': 1, b'TCGTATTTTT': 1, b'GGGGTGTCTT': 1, b'GNTGTCCTGA': 1, b'NTTTACAANT': 1, b'GGTACGTTAC': 1, b'TCTAAACNAN': 1, b'TATGGTNTTA': 1, b'GATTATCTTC': 1, b'CGANACCACC': 1, b'CTAACAACGG': 1, b'AGAGNCCTGA': 1, b'TTCTAATTTT': 1, b'TCTTNATCCT': 1, b'GAAANCATTA': 1, b'AAAATCGGNG': 1, b'CCCACCTATC': 1, b'TGAGGTGTCG': 1})}),
        'TTCATAAGGTAGGGAT':
        defaultdict(None, {
            'test3-CGTCGTAGCTGATCGTAGCTGAC': Counter({b'ATGACACCAT': 2, b'TTCTAATTAT': 2, b'TCCATGCAGT': 2, b'NATTTGGGGT': 2, b'AGAGGCCGTN': 2, b'ATANCGTAGA': 2, b'TGGCCCGTGA': 2, b'GAAATAATAG': 2, b'GGGTCGGCAC': 1, b'GTCCCAGTTC': 1, b'CTACCTTTGA': 1, b'TGCNCGTGAA': 1, b'GCAGNCACNG': 1, b'NATATCCCAT': 1, b'TCACANTATC': 1, b'TAAGAAGAAC': 1, b'GGCCCACGGT': 1, b'CACTAATCTT': 1, b'TNCCCAAGNT': 1, b'GCGATAACTN': 1, b'AGGGACTGGC': 1, b'CATTATGTAA': 1, b'CCATGGGACT': 1, b'ANNACCGGGT': 1, b'GACAACCCAA': 1, b'TCCAAATATT': 1, b'NTNANGGGNC': 1}),
            'test2-CGTACGTAGCCTAGC': Counter({b'GACAACCCAA': 2, b'CCAGTTTGTC': 2, b'CCCCGCGCAA': 2, b'TTCTAATTAT': 2, b'NGCANAGANN': 1, b'CCATAAACCA': 1, b'GTCAGTTGCG': 1, b'TTAGTAAAGT': 1, b'GTCCCAGTTC': 1, b'CGACNTGCCG': 1, b'AATATGNTTA': 1, b'AACCTTCGAG': 1, b'NATATCCCAT': 1, b'TCTAGGAAGT': 1, b'AACGNCGTTC': 1, b'GTNGGTGTCC': 1, b'TTAACTGTTA': 1, b'GGCTCTGCGN': 1, b'TACNGGGCGG': 1, b'AGGTAGNTCA': 1, b'GGTGCTATGG': 1, b'TNCCCAAGNT': 1, b'AGGGACTGGC': 1, b'AAGGTTCCAG': 1}),
            'test1-CGTAGCTCG': Counter({b'CGACNTGCCG': 2, b'AGAGGCCGTN': 2, b'TACNGGGCGG': 2, b'AACCTTCGAG': 2, b'TCCCCGCTCA': 2, b'AGGCCCTACT': 1, b'CATTATGTAA': 1, b'TTANCAGTCN': 1, b'ANTAGCNATG': 1, b'TCATGACCGA': 1, b'TAAACGCGGG': 1, b'ATANCGTAGA': 1, b'TCACANTATC': 1, b'AACGNCGTTC': 1, b'TTAACTGTTA': 1, b'TCATNACGGC': 1, b'TCCAAATATT': 1, b'ATCGANCTAC': 1, b'GGTGAACCAT': 1, b'CACTAATCTT': 1, b'CGATCGTAGA': 1, b'GGTANGACCT': 1, b'AGGGACTGGC': 1, b'ANTATTACCA': 1, b'CAACNTACAG': 1, b'CTCGCGTCTN': 1, b'CAGTAATCTT': 1, b'CACTTAGCAC': 1, b'GCTGTATCCG': 1, b'GGCCCACGGT': 1, b'GCAGNCACNG': 1, b'GCGATAACTN': 1})})}
    pytest.no_match = Counter()
    pytest.n = 200
    
    pytest.corrected_results = {
    'ACTGTTTTATTGGCCT':
    defaultdict(None, {
        'test3-CGTCGTAGCTGATCGTAGCTGAC': Counter({b'CATTAGTGGT': 2, b'ACACTCNTGT': 2, b'AAGCCTTTAN': 1, b'CCCACGGGTG': 1, b'AAACAAAAGA': 1, b'ACTTCTTAGC': 1, b'AGAGNCCTGA': 1, b'GACNATTGCA': 1, b'TTCAATNGGC': 1, b'CNCCTCGCTC': 1, b'AACGTCGTTT': 1, b'CNGNTATCTT': 1, b'GCTATTCTGG': 1, b'CCAGATTAGT': 1, b'GNTGTCCTGA': 1, b'TCTAAACNAN': 1, b'NAANNTTATA': 1, b'CCCCGGCAAT': 1, b'GATAGTNCGA': 1, b'TATGGTNTTA': 1, b'GTGGTGATAT': 1, b'GATTATCTTC': 1, b'TCGCTTGTAT': 1, b'CAGAACTTCT': 1, b'ACAGGCCACA': 1, b'TTACTTCCTC': 1, b'ATCTCTCAAA': 1, b'CACGTTANAG': 1, b'GCCCANATTA': 1, b'CGGGAGTGCG': 1}),
        'test1-CGTAGCTCG': Counter({b'CNCCTCGCTC': 2, b'CTGTGGGACA': 2, b'AAAATCGGNG': 2, b'AAGCCTTTAN': 2, b'TCTTNATCCT': 1, b'GTGGTGATAT': 1, b'TTGAATACTG': 1, b'ANTACATGTA': 1, b'GANCAGCCTC': 1, b'ATAGGGTTGC': 1, b'ATAACCCCTC': 1, b'CAGAACTTCT': 1, b'AGTATAATAC': 1, b'TATGCTTTGG': 1, b'TGAAATCTTG': 1, b'GGGTTAGCGC': 1, b'ACTTCTTAGC': 1, b'CTTTTNCCNG': 1, b'AATTCAAATC': 1, b'AGAGNCCTGA': 1, b'CATTAGTGGT': 1, b'GNTGTCCTGA': 1, b'GCTACAAGGG': 1, b'TAGTGCTGAG': 1, b'TCGGTCACAA': 1, b'AAACGCGCCG': 1, b'GAAANCATTA': 1, b'CCTACAAATA': 1, b'TTGCAGTGAT': 1, b'CGGGAGTGCG': 1, b'AAGTTTCTCT': 1, b'CGAGGGGCTG': 1}),
        'test2-CGTACGTAGCCTAGC': Counter({b'TCGTATGGCN': 2, b'CGGGAGTGCG': 2, b'CCAGGTGGGA': 2, b'TTACTTCCTC': 2, b'GANCAGCCTC': 1, b'CATTAGTGGT': 1, b'AATTCAAATC': 1, b'CACGTTANAG': 1, b'CNGNTATCTT': 1, b'NAANNTTATA': 1, b'TCANCCNCCG': 1, b'TCGTATTTTT': 1, b'GGGGTGTCTT': 1, b'GNTGTCCTGA': 1, b'NTTTACAANT': 1, b'GGTACGTTAC': 1, b'TCTAAACNAN': 1, b'TATGGTNTTA': 1, b'GATTATCTTC': 1, b'CGANACCACC': 1, b'CTAACAACGG': 1, b'AGAGNCCTGA': 1, b'TTCTAATTTT': 1, b'TCTTNATCCT': 1, b'GAAANCATTA': 1, b'AAAATCGGNG': 1, b'CCCACCTATC': 1, b'TGAGGTGTCG': 1})}),
    'TTCATAAGGTAGGGAT':
    defaultdict(None, {
        'test3-CGTCGTAGCTGATCGTAGCTGAC': Counter({b'ATGACACCAT': 2, b'TTCTAATTAT': 2, b'TCCATGCAGT': 2, b'NATTTGGGGT': 2, b'AGAGGCCGTN': 2, b'ATANCGTAGA': 2, b'TGGCCCGTGA': 2, b'GAAATAATAG': 2, b'GGGTCGGCAC': 1, b'GTCCCAGTTC': 1, b'CTACCTTTGA': 1, b'TGCNCGTGAA': 1, b'GCAGNCACNG': 1, b'NATATCCCAT': 1, b'TCACANTATC': 1, b'TAAGAAGAAC': 1, b'GGCCCACGGT': 1, b'CACTAATCTT': 1, b'TNCCCAAGNT': 1, b'GCGATAACTN': 1, b'AGGGACTGGC': 1, b'CATTATGTAA': 1, b'CCATGGGACT': 1, b'ANNACCGGGT': 1, b'GACAACCCAA': 1, b'TCCAAATATT': 1, b'NTNANGGGNC': 1}),
        'test2-CGTACGTAGCCTAGC': Counter({b'GACAACCCAA': 2, b'CCAGTTTGTC': 2, b'CCCCGCGCAA': 2, b'TTCTAATTAT': 2, b'NGCANAGANN': 1, b'CCATAAACCA': 1, b'GTCAGTTGCG': 1, b'TTAGTAAAGT': 1, b'GTCCCAGTTC': 1, b'CGACNTGCCG': 1, b'AATATGNTTA': 1, b'AACCTTCGAG': 1, b'NATATCCCAT': 1, b'TCTAGGAAGT': 1, b'AACGNCGTTC': 1, b'GTNGGTGTCC': 1, b'TTAACTGTTA': 1, b'GGCTCTGCGN': 1, b'TACNGGGCGG': 1, b'AGGTAGNTCA': 1, b'GGTGCTATGG': 1, b'TNCCCAAGNT': 1, b'AGGGACTGGC': 1, b'AAGGTTCCAG': 1}),
        'test1-CGTAGCTCG': Counter({b'CGACNTGCCG': 2, b'AGAGGCCGTN': 2, b'CACTAATCTT': 2, b'TACNGGGCGG': 2, b'AACCTTCGAG': 2, b'TCCCCGCTCA': 2, b'AGGCCCTACT': 1, b'CATTATGTAA': 1, b'TTANCAGTCN': 1, b'ANTAGCNATG': 1, b'TCATGACCGA': 1, b'TAAACGCGGG': 1, b'ATANCGTAGA': 1, b'TCACANTATC': 1, b'AACGNCGTTC': 1, b'TTAACTGTTA': 1, b'TCATNACGGC': 1, b'TCCAAATATT': 1, b'ATCGANCTAC': 1, b'GGTGAACCAT': 1, b'CGATCGTAGA': 1, b'GGTANGACCT': 1, b'AGGGACTGGC': 1, b'ANTATTACCA': 1, b'CAACNTACAG': 1, b'CTCGCGTCTN': 1, b'CACTTAGCAC': 1, b'GCTGTATCCG': 1, b'GGCCCACGGT': 1, b'GCAGNCACNG': 1, b'GCGATAACTN': 1})})}
    pytest.ordered_tags_map = OrderedDict({
        'test3-CGTCGTAGCTGATCGTAGCTGAC':0,
        'test2-CGTACGTAGCCTAGC':1,
        'test1-CGTAGCTCG': 3,
        'unmapped': 4
        })

@pytest.mark.dependency()
def test_classify_reads_multi_process(data):
    (results, no_match, n) = processing.map_reads(
        pytest.correct_R1_path,
        pytest.correct_R2_path,
        pytest.chunk_size,
        pytest.tags,
        pytest.barcode_slice,
        pytest.umi_slice,
        pytest.first_line,
        pytest.correct_whitelist,
        pytest.legacy,
        pytest.debug,
        pytest.start_trim,
        pytest.maximum_distance)
    for cell_barcode in results.keys():
        print(cell_barcode)
        for TAG in results[cell_barcode]:
            for UMI in results[cell_barcode][TAG]:
                print(results[cell_barcode][TAG][UMI])
                print(pytest.results[cell_barcode][TAG][UMI])
                #assert results[cell_barcode][TAG][UMI] == pytest.results[cell_barcode][TAG][UMI]
    #== (pytest.results, pytest.no_match, pytest.n)

@pytest.mark.dependency(depends=['test_classify_reads_multi_process'])
def test_correct_umis(data):
    temp = processing.correct_umis(pytest.results, 2)
    results = temp[0]
    n_corrected = temp[1]
    for cell_barcode in results.keys():
        for TAG in results[cell_barcode]:
            assert len(results[cell_barcode][TAG]) == len(pytest.corrected_results[cell_barcode][TAG])
            assert sum(results[cell_barcode][TAG].values()) == sum(pytest.corrected_results[cell_barcode][TAG].values())
    assert n_corrected == 1

@pytest.mark.dependency(depends=['test_correct_umis'])
def test_generate_sparse_matrices(data):
    (umi_results_matrix, read_results_matrix) = processing.generate_sparse_matrices(
        pytest.corrected_results, pytest.ordered_tags_map,
        set(['ACTGTTTTATTGGCCT','TTCATAAGGTAGGGAT'])
        )
    assert umi_results_matrix.shape == (4,2)
    assert read_results_matrix.shape == (4,2)
    read_results_matrix = read_results_matrix.tocsr()
    total_reads = 0
    for i in range(read_results_matrix.shape[0]):
        for j in range(read_results_matrix.shape[1]):
            total_reads += read_results_matrix[i,j]
    assert total_reads == 200