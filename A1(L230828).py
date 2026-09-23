def pattern_Count(text, pattern):
    count = 0
    pattern_length = len(pattern)
    for i in range(len(text) - pattern_length + 1):
        if text.startswith(pattern, i):
            count += 1
    return count

print(pattern_Count("ACGTACGTACGTGT", "GT"))