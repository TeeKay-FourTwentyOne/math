# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Purpose

This is an experimental project for using AI systems to attack open math problems, starting with Ramsey numbers for book graphs. The focus is on tooling and infrastructure to support mathematical research.

## Current Problem: R(B_{n-1}, B_n) = 4n - 1

**Book graph B_n**: K_2 + K̄_n (n triangles sharing a common edge)

**Goal**: Prove R(B_{n-1}, B_n) = 4n - 1 for all n

**Known results**:
- Upper bound ≤ 4n - 1 proven (Rousseau & Sheehan, 1978)
- Lower bound verified for n ≤ 21 and when 2n - 1 is a prime power ≡ 1 (mod 4)
- Constructions use 2-block circulant graphs and Paley graph generalizations

**Key references**: See ramsey-book-graphs.pdf for full problem statement and citations.

## Project Structure

- `ramsey-book-graphs.pdf` - Problem statement from FrontierMath/Epoch AI
- `*.csv` files - Supplementary data (ignore for now)
