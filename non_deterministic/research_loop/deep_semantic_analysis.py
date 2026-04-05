"""
Deep semantic analysis of word-vs-word competing tokens.

word_vs_word is the largest category (39.7%) of fragile position competitions.
Let's understand what kind of word substitutions are happening.
"""

import torch
import numpy as np
import json
import os
from collections import Counter, defaultdict

from transformers import AutoTokenizer
tokenizer = AutoTokenizer.from_pretrained("EleutherAI/pythia-410m")

LOGIT_DIR = "C:/Projects/math/non_deterministic/experiment_0/logit_traces"
PROMPTS_FILE = "C:/Projects/math/non_deterministic/prompts/candidates_2000.jsonl"
MARGIN_THRESHOLD = 0.1

prompts = {}
with open(PROMPTS_FILE) as f:
    for line in f:
        d = json.loads(line)
        prompts[d['idx']] = d

def decode(tid):
    return tokenizer.decode([tid])

def classify_token(text):
    text_stripped = text.strip()
    if not text_stripped:
        return "whitespace"
    if text_stripped.isdigit():
        return "digit"
    if all(c in '.,;:!?-\'"()[]{}' for c in text_stripped):
        return "punctuation"
    if text_stripped.startswith('\n'):
        return "newline"
    if text_stripped.isalpha():
        if text_stripped[0].isupper():
            return "word_capitalized"
        return "word_lower"
    if any(c.isdigit() for c in text_stripped) and any(c.isalpha() for c in text_stripped):
        return "alphanumeric"
    return "other"

# Collect all word-vs-word pairs
word_pairs = []  # (prompt_idx, pos, tok1_text, tok2_text, margin, category)

for i in range(1999):
    path = os.path.join(LOGIT_DIR, f"prompt_{i:04d}.pt")
    if not os.path.exists(path):
        continue
    data = torch.load(path, map_location='cpu', weights_only=False)
    top5_vals = data['top5_values']
    top5_idx = data['top5_indices']
    margins = top5_vals[:, 0] - top5_vals[:, 1]

    for pos in range(len(margins)):
        if margins[pos] < MARGIN_THRESHOLD:
            t1 = decode(top5_idx[pos, 0].item())
            t2 = decode(top5_idx[pos, 1].item())
            c1 = classify_token(t1)
            c2 = classify_token(t2)
            if "word" in c1 and "word" in c2:
                t1s = t1.strip().lower()
                t2s = t2.strip().lower()
                if t1s != t2s:  # not case variant
                    word_pairs.append((i, pos, t1, t2, margins[pos].item(), prompts[i]['category']))

print(f"Total word-vs-word fragile positions: {len(word_pairs)}")

# Classify the semantic relationship more finely
def fine_classify_word_pair(t1, t2):
    """Fine-grained classification of word pairs."""
    w1 = t1.strip().lower()
    w2 = t2.strip().lower()

    # Prefix/suffix relationship
    if w1.startswith(w2) or w2.startswith(w1):
        return "prefix_substring"

    # Article/determiner swaps
    articles = {'a', 'an', 'the', 'this', 'that', 'these', 'those', 'some', 'any', 'no', 'each', 'every'}
    if w1 in articles and w2 in articles:
        return "article_swap"

    # Pronoun swaps
    pronouns = {'i', 'me', 'my', 'mine', 'we', 'us', 'our', 'ours',
                'you', 'your', 'yours', 'he', 'him', 'his', 'she', 'her', 'hers',
                'it', 'its', 'they', 'them', 'their', 'theirs'}
    if w1 in pronouns and w2 in pronouns:
        return "pronoun_swap"

    # Preposition swaps
    preps = {'in', 'on', 'at', 'to', 'for', 'with', 'by', 'from', 'of', 'about',
             'into', 'through', 'during', 'before', 'after', 'above', 'below',
             'between', 'under', 'over', 'among', 'around', 'within'}
    if w1 in preps and w2 in preps:
        return "preposition_swap"

    # Conjunction/connector swaps
    conj = {'and', 'or', 'but', 'nor', 'yet', 'so', 'however', 'although', 'while', 'whereas'}
    if w1 in conj and w2 in conj:
        return "conjunction_swap"

    # Function word vs content word
    function_words = articles | pronouns | preps | conj | {'is', 'was', 'are', 'were', 'be', 'been', 'being',
        'has', 'have', 'had', 'do', 'does', 'did', 'will', 'would', 'shall', 'should',
        'can', 'could', 'may', 'might', 'must', 'not', 'also', 'just', 'very', 'only',
        'quite', 'rather', 'still', 'even', 'more', 'most', 'less', 'as', 'than',
        'if', 'then', 'when', 'where', 'how', 'what', 'which', 'who', 'whom', 'whose'}

    both_function = w1 in function_words and w2 in function_words
    both_content = w1 not in function_words and w2 not in function_words

    # Tense/form variants
    # Check if one is a common morphological variant
    if (w1 + 's' == w2 or w2 + 's' == w1 or
        w1 + 'ed' == w2 or w2 + 'ed' == w1 or
        w1 + 'ing' == w2 or w2 + 'ing' == w1 or
        w1 + 'er' == w2 or w2 + 'er' == w1 or
        w1 + 'ly' == w2 or w2 + 'ly' == w1 or
        w1 + 'es' == w2 or w2 + 'es' == w1):
        return "morphological_variant"

    # Tense swaps (common verb pairs)
    tense_pairs = {
        ('is', 'was'), ('are', 'were'), ('has', 'had'), ('have', 'had'),
        ('do', 'did'), ('does', 'did'), ('will', 'would'), ('can', 'could'),
        ('shall', 'should'), ('may', 'might'), ('is', 'are'), ('was', 'were'),
        ('has', 'have'), ('does', 'do'),
    }
    if (w1, w2) in tense_pairs or (w2, w1) in tense_pairs:
        return "tense_swap"

    # Sentence starters (common first words of sentences)
    starters = {'the', 'a', 'it', 'this', 'there', 'he', 'she', 'they', 'we',
                'in', 'on', 'for', 'however', 'although', 'but', 'so',
                'one', 'some', 'many', 'most', 'all', 'each', 'every'}

    if both_function:
        return "function_word_swap"
    if both_content:
        return "content_word_swap"

    return "mixed_function_content"

fine_cats = Counter()
fine_examples = defaultdict(list)

for prompt_idx, pos, t1, t2, margin, cat in word_pairs:
    fcat = fine_classify_word_pair(t1, t2)
    fine_cats[fcat] += 1
    if len(fine_examples[fcat]) < 8:
        fine_examples[fcat].append((prompt_idx, pos, t1, t2, margin, cat))

print("\n" + "=" * 70)
print("Fine-grained word pair classification:")
print("=" * 70)
for fcat, count in fine_cats.most_common():
    pct = 100 * count / len(word_pairs)
    print(f"\n{fcat}: {count} ({pct:.1f}%)")
    for ex in fine_examples[fcat][:5]:
        prompt_idx, pos, t1, t2, margin, cat = ex
        print(f"  [{cat}] p{prompt_idx} @{pos}: {repr(t1):15s} vs {repr(t2):15s} (m={margin:.4f})")

# Now the critical question: among content_word_swap, are these semantically
# related (synonyms, antonyms) or completely different words?
print("\n" + "=" * 70)
print("Content word swaps: detailed examples (first 50)")
print("=" * 70)

content_swaps = [(pi, pos, t1, t2, m, cat) for pi, pos, t1, t2, m, cat in word_pairs
                 if fine_classify_word_pair(t1, t2) == 'content_word_swap']

# Manual categorization of the first 50 content word swaps
print(f"\nTotal content word swaps: {len(content_swaps)}")
for i, (prompt_idx, pos, t1, t2, margin, cat) in enumerate(content_swaps[:50]):
    # Get some context
    path = os.path.join(LOGIT_DIR, f"prompt_{prompt_idx:04d}.pt")
    data = torch.load(path, map_location='cpu', weights_only=False)
    gen_ids = data['generated_ids']

    # Decode context around the fragile position
    input_len = len(data['input_ids'])
    gen_start = pos  # position in generation (0-indexed)
    # The generated token at this position
    actual_pos_in_seq = input_len + gen_start
    # Get context: 5 tokens before
    ctx_start = max(0, actual_pos_in_seq - 5)
    ctx_end = min(len(gen_ids), actual_pos_in_seq + 3)
    context = tokenizer.decode(gen_ids[ctx_start:actual_pos_in_seq])
    after = tokenizer.decode(gen_ids[actual_pos_in_seq:ctx_end])

    print(f"\n{i+1}. [{cat}] p{prompt_idx} @pos{pos}: {repr(t1)} vs {repr(t2)} (m={margin:.4f})")
    print(f"   Prompt: '{prompts[prompt_idx]['text'][:80]}'")
    print(f"   Context: '...{context}[{t1.strip()}|{t2.strip()}]{after}...'")


# Now analyze: at early positions, the word swaps are often about the discourse
# direction. Let's see what fraction of content swaps at pos 0-10 vs later.
print("\n" + "=" * 70)
print("Content word swap position distribution")
print("=" * 70)
pos_bins = [(0, 7), (8, 15), (16, 31), (32, 63), (64, 127), (128, 255)]
for lo, hi in pos_bins:
    count = sum(1 for _, p, _, _, _, _ in content_swaps if lo <= p <= hi)
    pct = 100 * count / len(content_swaps) if content_swaps else 0
    print(f"  pos {lo:3d}-{hi:3d}: {count:4d} ({pct:.1f}%)")


# Finally: what about the `' was' vs ' is'` pattern?
# This was seen at position 16 and could be tense confusion
print("\n" + "=" * 70)
print("Tense confusion examples ('was' vs 'is' and similar):")
print("=" * 70)
tense_count = 0
for prompt_idx, pos, t1, t2, margin, cat in word_pairs:
    w1 = t1.strip().lower()
    w2 = t2.strip().lower()
    if (w1, w2) in [('was', 'is'), ('is', 'was'), ('are', 'were'), ('were', 'are'),
                      ('has', 'had'), ('had', 'has'), ('have', 'had'), ('had', 'have'),
                      ('does', 'did'), ('did', 'does'), ('do', 'did'), ('did', 'do')]:
        tense_count += 1
        if tense_count <= 10:
            print(f"  [{cat}] p{prompt_idx} @{pos}: {repr(t1):12s} vs {repr(t2):12s} (m={margin:.4f})")
            print(f"    '{prompts[prompt_idx]['text'][:80]}'")
print(f"\nTotal tense confusion pairs: {tense_count} ({100*tense_count/len(word_pairs):.1f}% of word pairs)")
