package edu.ukma.bioinfcalcweb.service;

import edu.ukma.bioinfcalcweb.util.IRandomProvider;
import org.springframework.stereotype.Service;

import java.util.*;

@Service
public class MotifSrv {
    private final IRandomProvider randomProvider;

    public MotifSrv(IRandomProvider randomProvider) {
        this.randomProvider = randomProvider;
    }

    public static List<String> greedyMotifSearch(String[] dna, int k) {
        List<String> bestMotifs = new ArrayList<>();
        for (String strand : dna)
            bestMotifs.add(strand.substring(0, k));
        int n = dna[0].length();
        for (int i = 0; i < n - k + 1; i++) {
            List<String> motifs = new ArrayList<>();
            motifs.add(dna[0].substring(i, i + k));
            for (int j = 1; j < dna.length; j++) {
                Map<Character, List<Double>> profile = profileWithPseudocounts(motifs.subList(0, j));
                motifs.add(profileMostProbableKmer(dna[j], k, profile));
            }
            if (score(motifs) < score(bestMotifs)) bestMotifs = motifs;
        }
        return bestMotifs;
    }

    public List<String> randomizedMotifSearch(String[] dna, int k) {
        List<String> bestMotifs = randomMotifs(dna, k);
        while (true) {
            Map<Character, List<Double>> profile = profileWithPseudocounts(bestMotifs);
            List<String> m = findMotifs(profile, dna);
            if (score(m) > score(bestMotifs))
                bestMotifs = m;
            else return bestMotifs;
        }
    }

    private List<String> randomMotifs(String[] dna, int k) {
        List<String> motifs = new ArrayList<>();
        for (String segment : dna) {
            int start = randomProvider.nextInt(dna[0].length() - k + 1);
            motifs.add(segment.substring(start, start + k));
        }
        int score = dna.length * k;
        while (true) {
            Map<Character, List<Double>> profile = profileWithPseudocounts(motifs);
            motifs = findMotifs(profile, dna);
            int newScore = score(motifs);
            if (newScore < score) score = newScore;
            else break;
        }
        return motifs;
    }

    private static List<String> findMotifs(Map<Character, List<Double>> profile, String[] dna) {
        List<String> motifs = new ArrayList<>();
        for (String strand : dna)
            motifs.add(profileMostProbableKmer(strand, profile.get('A').size(), profile));
        return motifs;
    }

    public List<String> gibbsSamplerSearch(String[] dna, int k, int n) {
        List<String> bestMotifs = randomMotifs(dna, k);
        for (int i = 0; i < n; i++) {
            int searchStrIdx = randomProvider.nextInt(dna.length);
            List<String> investedMotifs = new ArrayList<>(Arrays.asList(dna));
            investedMotifs.remove(searchStrIdx);
            Map<Character, List<Double>> profile = profileWithPseudocounts(investedMotifs);
            String kmer = profileGeneratedStringDna(dna[searchStrIdx], profile, k);
            List<String> updMotifs = new ArrayList<>(bestMotifs);
            updMotifs.set(searchStrIdx, kmer);
            if (score(updMotifs) < score(bestMotifs)) bestMotifs = updMotifs;
        }
        return bestMotifs;
    }

    private String profileGeneratedStringDna(String text, Map<Character, List<Double>> profile, int k) {
        Map<String, Double> probabilities = new HashMap<>();
        for (int i = 0; i < text.length() - k + 1; i++) {
            String dnaStr = text.substring(i, i + k);
            probabilities.put(dnaStr, pr(dnaStr, profile));
        }
        return weightedDie(normalize(probabilities));
    }

    public String weightedDie(Map<String, Double> prob) {
        double selectedProb = randomProvider.nextDouble();
        double accum = 0;
        for (String kmer : prob.keySet()) {
            accum += prob.get(kmer);
            if (selectedProb <= accum) return kmer;
        }
        // probabilities may add up to 0.998, but random chooses 0.999
        return prob.keySet().toArray(new String[0])[prob.size() - 1];
    }

    private static Map<String, Double> normalize(Map<String, Double> probabilities) {
        double coef = 1.0 / probabilities.values().stream().reduce(0.0, Double::sum);
        Map<String, Double> normProfile = new HashMap<>(probabilities.size());
        for (Map.Entry<String, Double> entry : probabilities.entrySet())
            normProfile.put(entry.getKey(), entry.getValue() * coef);
        return normProfile;
    }

    private static Map<Character, List<Double>> profileWithPseudocounts(List<String> motifs) {
        Map<Character, List<Double>> countMat = new HashMap<>();
        int k = motifs.get(0).length();
        for (char symbol : "ACGT".toCharArray()) {
            countMat.put(symbol, new ArrayList<>());
            for (int i = 0; i < k; i++)
                countMat.get(symbol).add(1.0);
        }
        for (String motif : motifs)
            for (int j = 0; j < k; j++) {
                char symbol = motif.charAt(j);
                countMat.get(symbol).set(j, countMat.get(symbol).get(j) + 1);
            }
        for (char nucKey : countMat.keySet())
            for (int i = 0; i < motifs.get(0).length(); i++)
                countMat.get(nucKey).set(i, countMat.get(nucKey).get(i) / (motifs.size() + 4));
        return countMat;
    }

    private static int score(List<String> motifs) {
        String cons = consensus(motifs);
        int score = 0;
        for (String motif : motifs)
            for (int j = 0; j < motifs.get(0).length(); j++)
                if (motif.charAt(j) != cons.charAt(j))
                    ++score;
        return score;
    }

    private static String profileMostProbableKmer(String text, int k, Map<Character, List<Double>> profile) {
        double maxProb = 0;
        Integer resStartIdx = null;
        for (int i = 0; i < text.length() - k + 1; i++) {
            double nextProb = pr(text.substring(i, i + k), profile);
            if (nextProb > maxProb) {
                maxProb = nextProb;
                resStartIdx = i;
            }
        }
        return resStartIdx != null ? text.substring(resStartIdx, resStartIdx + k) : text.substring(0, k);
    }

    private static double pr(String text, Map<Character, List<Double>> profile) {
        double init = 1;
        for (int i = 0; i < text.length(); i++)
            init *= profile.get(text.charAt(i)).get(i);
        return init;
    }

    private static String consensus(List<String> motifs) {
        Map<Character, List<Integer>> count = countNucsPositions(motifs);
        StringBuilder consensus = new StringBuilder();
        for (int i = 0; i < motifs.get(0).length(); i++) {
            int m = 0;
            char frequentSymbol = 0;
            for (char symbol : "ACTG".toCharArray()) {
                if (count.get(symbol).get(i) > m) {
                    m = count.get(symbol).get(i);
                    frequentSymbol = symbol;
                }
            }
            consensus.append(frequentSymbol);
        }
        return consensus.toString();
    }

    private static Map<Character, List<Integer>> countNucsPositions(List<String> motifs) {
        int k = motifs.get(0).length();
        Map<Character, List<Integer>> count = new HashMap<>();
        for (char symbol : "ACGT".toCharArray()) {
            count.put(symbol, new ArrayList<>());
            for (int i = 0; i < k; i++)
                count.get(symbol).add(0);
        }
        for (String motif : motifs)
            for (int j = 0; j < k; j++) {
                char symbol = motif.charAt(j);
                count.get(symbol).set(j, count.get(symbol).get(j) + 1);
            }
        return count;
    }

}
