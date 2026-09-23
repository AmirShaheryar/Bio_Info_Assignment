# Q1

def pattern_Count(text, pattern):
    count = 0
    pattern_length = len(pattern)
    for i in range(len(text) - pattern_length + 1):
        if text.startswith(pattern, i):
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

print(frequent_words("CGTTTTGAACATTTTCAACAAGTTTTGCAACATTTT ", 4))