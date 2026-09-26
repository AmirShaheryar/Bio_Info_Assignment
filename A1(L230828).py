from itertools import product

# Q1

def pattern_Count(text, pattern):
    count = 0

    for i in range(len(text) - len(pattern) + 1):
        if text[i:i+len(pattern)] == pattern:
            count += 1

    return count

print(pattern_Count("ACGTACGTACGTGT", "GT"))

#Q2

def frequent_words(text, k):
    freq_map = {}
    for i in range(len(text) - k + 1):
        kmer = text[i:i+k]
        if kmer in freq_map:
            freq_map[kmer] += 1
        else:
            freq_map[kmer] = 1


    max_count = max(freq_map.values())
    
    frequent_kmers = []
    for kmer, count in freq_map.items():
        if count == max_count:
            frequent_kmers.append(kmer)
            
    return frequent_kmers

print(frequent_words("CGTTTTGAACATTTTCAACAAGTTTTGCAACATTTTAACA ", 4))

#Q3

def Minimum_Skew(text):
    skew = 0
    min_skew = 0
    min_positions = [0]
    
    for i in range(len(text)):
        if text[i] == 'G':
            skew += 1
        elif text[i] == 'C':
            skew -= 1
        
        if skew < min_skew:
            min_skew = skew
            min_positions = [i + 1]  
        elif skew == min_skew:
            min_positions.append(i + 1)  
    
    return min_positions

print(Minimum_Skew("CCGGCCGG"))
print("\n")
print(Minimum_Skew("CCGCCGG"))

#Q4

def frequent_words_with_mismatches(text, k, d):
    from collections import defaultdict

    def neighbors(pattern, d):
        if d == 0:
            return {pattern}
        if len(pattern) == 0:
            return {""}
        
        neighborhood = set()
        suffix_neighbors = neighbors(pattern[1:], d)
        for text in suffix_neighbors:
            if hamming_distance(pattern[1:], text) < d:
                for nucleotide in "ACGT":
                    neighborhood.add(nucleotide + text)
            else:
                neighborhood.add(pattern[0] + text)
        return neighborhood

    def hamming_distance(s1, s2):
        return sum(el1 != el2 for el1, el2 in zip(s1, s2))

    freq_map = defaultdict(int)
    for i in range(len(text) - k + 1):
        kmer = text[i:i+k]
        neighborhood = neighbors(kmer, d)
        for neighbor in neighborhood:
            freq_map[neighbor] += 1

    max_count = max(freq_map.values())
    frequent_kmers = [kmer for kmer, count in freq_map.items() if count == max_count]

    return frequent_kmers

print(frequent_words_with_mismatches("ACGTTGCATGTCGCATGATGCATGAGAGCT", 4, 1))

#Q5

from itertools import product

def median_string(dna, k):

    def hamming_distance(s1, s2):
        count = 0

        for i in range(len(s1)):
            if s1[i] != s2[i]:
                count += 1

        return count

    def all_kmers(k):
        kmers = []

        for p in product("ACGT", repeat=k):
            kmers.append("".join(p))

        return kmers

    kmers = all_kmers(k)

    min_distance = float("inf")
    median_kmer = None

    for kmer in kmers:

        total_distance = 0

        for seq in dna.split():

            min_seq_distance = float("inf")

            for i in range(len(seq) - k + 1):

                substring = seq[i:i+k]

                distance = hamming_distance(kmer, substring)

                if distance < min_seq_distance:
                    min_seq_distance = distance

            total_distance += min_seq_distance

        if total_distance < min_distance:
            min_distance = total_distance
            median_kmer = kmer

    return median_kmer

print(median_string("AAATTGACGCAT GACGACCACGTT CGTCAGCGCCTG GCTGAGCACCGG AGTACGGGACAG", 3))

#Q6

def greedy_motif_search(dna, k, t):

    if isinstance(dna, str):
        dna = dna.split()

    def profile_most_probable_kmer(text, k, profile):

        best_kmer = ""
        best_probability = -1

        for i in range(len(text) - k + 1):

            kmer = text[i:i+k]
            probability = 1

            for pos in range(k):
                base = kmer[pos]
                probability = probability * profile[base][pos]

            if probability > best_probability:
                best_probability = probability
                best_kmer = kmer

        return best_kmer

    def profile_with_pseudocounts(motifs):

        profile = {
            "A": [1] * k,
            "C": [1] * k,
            "G": [1] * k,
            "T": [1] * k
        }

        for motif in motifs:

            for i in range(k):

                base = motif[i]
                profile[base][i] += 1

        for base in "ACGT":

            total = sum(profile[base])

            for i in range(k):
                profile[base][i] = profile[base][i] / total

        return profile

    def score(motifs):

        total = 0

        for column in zip(*motifs):

            max_count = 0

            for base in "ACGT":

                count = column.count(base)

                if count > max_count:
                    max_count = count

            total = total + (len(column) - max_count)

        return total

    best_motifs = []

    for i in range(len(dna[0]) - k + 1):
        best_motifs.append(dna[0][i:i+k])

    best_motifs = best_motifs[:t]

    for start in range(len(dna[0]) - k + 1):

        motif = dna[0][start:start+k]

        motifs = [motif]

        for i in range(1, t):

            profile = profile_with_pseudocounts(motifs)

            best = profile_most_probable_kmer(
                dna[i], k, profile
            )

            motifs.append(best)

        if score(motifs) < score(best_motifs):
            best_motifs = motifs

    return best_motifs

print(greedy_motif_search("GGCGTTCAGGCA AAGAATCAGTCA CAAGGAGTTCGC CACGTCAATCAC CAATAATATTCG", 3, 5))

#Q7

def overlap_graph(reads):
    graph = {}
    for read in reads:
        suffix = read[1:]
        for other_read in reads:
            if read != other_read and other_read.startswith(suffix):
                if read not in graph:
                    graph[read] = []
                graph[read].append(other_read)
    return graph

print(overlap_graph(["AAG", "AGA", "ATT", "CTA", "CTC", "GAT", "TAC", "TCT", "TCT", "TTC"]))

#Q8

def de_bruijn_kmers(kmers):
    graph = {}
    for kmer in kmers:
        prefix = kmer[:-1]
        suffix = kmer[1:]
        if prefix not in graph:
            graph[prefix] = []
        graph[prefix].append(suffix)
    return graph

print(de_bruijn_kmers(["GAGG", "CAGG", "GGGG", "GGGA", "CAGG", "AGGG", "GGAG"]))

#Q9

def stringReconstructionReadPairs(pairs, k, d):
    from collections import defaultdict

    def build_de_bruijn_graph(pairs):
        graph = defaultdict(list)
        for a, b in pairs:
            prefix = (a[:-1], b[:-1])
            suffix = (a[1:], b[1:])
            graph[prefix].append(suffix)
        return graph

    def find_eulerian_path(graph):
        from collections import defaultdict, deque

        in_degree = defaultdict(int)
        out_degree = defaultdict(int)

        for node in graph:
            out_degree[node] += len(graph[node])
            for neighbor in graph[node]:
                in_degree[neighbor] += 1

        start_node = None
        end_node = None
        for node in set(in_degree.keys()).union(set(out_degree.keys())):
            if out_degree[node] - in_degree[node] == 1:
                start_node = node
            elif in_degree[node] - out_degree[node] == 1:
                end_node = node

        if start_node is None:
            start_node = next(iter(graph))

        stack = [start_node]
        path = []

        while stack:
            current_node = stack[-1]
            if current_node in graph and graph[current_node]:
                next_node = graph[current_node].pop()
                stack.append(next_node)
            else:
                path.append(stack.pop())

        return path[::-1]

    def reconstruct_string(path, k, d):
        prefix_string = path[0][0]
        suffix_string = path[0][1]

        for i in range(1, len(path)):
            prefix_string += path[i][0][-1]
            suffix_string += path[i][1][-1]

        return prefix_string + suffix_string[-(k + d):]

    de_bruijn_graph = build_de_bruijn_graph(pairs)
    eulerian_path = find_eulerian_path(de_bruijn_graph)
    reconstructed_string = reconstruct_string(eulerian_path, k, d)

    return reconstructed_string

print(stringReconstructionReadPairs([("ACAC", "CTCT"), ("ACAT", "CTCA"), ("CACA", "TCTC"), ("GACA", "TCTC")], 4, 2))