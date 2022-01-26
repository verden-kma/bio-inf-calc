package edu.ukma.bioinfcalcweb.service;

import com.google.common.collect.*;
import com.google.common.primitives.Ints;
import edu.ukma.bioinfcalcweb.util.Constants;
import org.springframework.data.util.Pair;

import java.util.*;
import java.util.stream.Collectors;

public class PeptideReconstructionSrv {

    public static Collection<int[]> perfectSpectrumPeptides(final int[] spectrum) {
        Map<Integer, Integer> spectrumSet = new HashMap<>(); // maybe use Multiset
        for (int i = 1; i < spectrum.length; i++) {
            spectrumSet.compute(spectrum[i], (acid, num) -> num == null ? 1 : num + 1);
        }
        Set<Integer> basicAcids = spectrumSet.keySet();
        basicAcids.retainAll(Constants.AMINO_ACID_MASS.values());

        BiMap<Integer, Integer> basicOrder = HashBiMap.create(basicAcids.size());
        for (Integer acidMass : basicAcids)
            basicOrder.put(basicOrder.size(), acidMass);

        // Queue is a collection of consistent peptide abstractions
        // fst pair = a peptide denoted by an array representing masses of its amino acids
        // snd pair = num of occurrences of each basic consistent amino acid
        Queue<Pair<int[], int[]>> consistentPeptides = new ArrayDeque<>();
        List<int[]> peptides = new ArrayList<>();
        basicAcids.forEach(acid -> {
            int[] freqCache = new int[basicAcids.size()];
            freqCache[basicOrder.inverse().get(acid)] = 1;
            consistentPeptides.add(Pair.of(new int[]{acid}, freqCache));
        });
        while (!consistentPeptides.isEmpty()) {
            Pair<int[], int[]> currPeptide = consistentPeptides.poll();
            for (Integer acidMass : basicAcids) {
                int[] acidsCache = currPeptide.getSecond();
                int alreadyPresent = acidsCache[basicOrder.inverse().get(acidMass)];
                if (spectrumSet.get(acidMass) == alreadyPresent) continue;

                int[] nextPeptide = Arrays.copyOf(currPeptide.getFirst(), currPeptide.getFirst().length + 1);
                nextPeptide[nextPeptide.length - 1] = acidMass;

                if (spectrumSet.get(acidMass) == alreadyPresent + 1) {
                    int peptideMass = 0;
                    for (int i = 0; i < acidsCache.length; i++) {
                        peptideMass += acidsCache[i] * basicOrder.get(i);
                    }
                    peptideMass += acidMass;
                    if (peptideMass == spectrum[spectrum.length - 1]
                            && Arrays.equals(spectrum, theoreticalPeptideSpectrum(nextPeptide, true))) {
                        peptides.add(nextPeptide);
                        continue;
                    }
                }
                int[] cacheUpd = Arrays.copyOf(acidsCache, acidsCache.length);
                cacheUpd[basicOrder.inverse().get(acidMass)]++;
                consistentPeptides.add(Pair.of(nextPeptide, cacheUpd));
            }
        }
        return peptides;
    }

    public static int[] theoreticalPeptideSpectrum(final String peptide, boolean isCircular) {
        SortedMultiset<Integer> spectrum = TreeMultiset.create();
        spectrum.add(0);
        int accum = 0;
        for (int i = 0; i < peptide.length(); i++) {
            accum = 0;
            for (int j = i; j < peptide.length() + (isCircular ? i - 1 : 0); j++) {
                accum += Constants.AMINO_ACID_MASS.get(peptide.charAt(j - (isCircular && j >= peptide.length() ? peptide.length() : 0)));
                spectrum.add(accum);
            }
        }
        if (peptide.length() > 1)
            spectrum.add(accum + Constants.AMINO_ACID_MASS.get(peptide.charAt(peptide.length() - 2)));
        return Ints.toArray(spectrum);
    }

    public static int[] theoreticalPeptideSpectrum(final int[] peptide, boolean isCircular) {
        SortedMultiset<Integer> spectrum = TreeMultiset.create();
        spectrum.add(0);
        int accum = 0;
        for (int i = 0; i < peptide.length; i++) {
            accum = 0;
            for (int j = i; j < peptide.length + (isCircular ? i - 1 : 0); j++) {
                accum += peptide[j - (isCircular && j >= peptide.length ? peptide.length : 0)];
                spectrum.add(accum);
            }
        }
        if (peptide.length > 1)
            spectrum.add(accum + peptide[peptide.length - 2]);
        return Ints.toArray(spectrum);
    }

    public static int doScorePeptide(final int[] pSpec, final int[] spectrum, boolean isCircular) {
        int score = 0;
        for (int i = 0, j = 0; i < pSpec.length && j < spectrum.length; ) {
            if (pSpec[i] == spectrum[j]) {
                i++;
                j++;
                score++;
            } else if (pSpec[i] > spectrum[j]) j++;
            else i++;
        }
        return score;
    }

    public static int scorePeptide(final String peptide, final int[] spectrum, boolean isCircular) {
        return doScorePeptide(theoreticalPeptideSpectrum(peptide, isCircular), spectrum, isCircular);
    }

    public static int scorePeptide(final int[] peptide, final int[] spectrum, boolean isCircular) {
        return doScorePeptide(theoreticalPeptideSpectrum(peptide, isCircular), spectrum, isCircular);
    }

    public static List<int[]> leaderboardCyclopeptideSequencing(final int[] spectrum, final int[] acidMasses, final int n) {
        final Multiset<Integer> spectrumIdx = HashMultiset.create();
        for (int mass : spectrum) spectrumIdx.add(mass);

        ArrayDeque<PeptideCandidate> candidates = new ArrayDeque<>();
        candidates.add(new PeptideCandidate());
        List<Leader> leaders = new ArrayList<>();
        TreeMap<Integer, Integer> scoreLen = new TreeMap<>();
        scoreLen.put(0, Integer.MAX_VALUE);

        while (!candidates.isEmpty()) {
            PeptideCandidate candidate = candidates.poll();
            for (Integer acidMass : acidMasses) {
                if (candidate.mass + acidMass > spectrum[spectrum.length - 1]) continue;
                int scoreInc = 0;
                // duplicate entries?
                // no - because next entry is a superset of a previous, so it is guaranteed to have higher mass
                int[] outgrowth = candidate.deltaSpectrum(acidMass);
                for (int spMass : outgrowth) {
                    if (spectrumIdx.count(spMass) - candidate.localSpectrum.count(spMass) > 0) scoreInc++;
                }
                int offeredScore = candidate.score + scoreInc;
                /*
                 * ignore candidates with the score lower, the such of the worst leader's,
                 * but let them through during initialization phase (in case noise at the beginning obstruct potentially
                 * optimal peptides' beginnings)
                 * also ignore candidates that offer the same score as their shorter counterparts
                 * IF ABCDE is correct and ABCm is already present, ignore ABCmp etc.
                 * */
                if (offeredScore < scoreLen.firstKey() && scoreLen.size() == n ||
                        scoreLen.containsKey(offeredScore) && scoreLen.get(offeredScore) < candidate.acids.length + 1)
                    continue;
                PeptideCandidate nextCandidate = candidate.ascend(acidMass, outgrowth, scoreInc);

                if (nextCandidate.mass == spectrum[spectrum.length - 1]) {
                    Multiset<Integer> cyclicPart = nextCandidate.cyclicSpectrumAddition();
                    int cyclicScore = nextCandidate.score;
                    for (Multiset.Entry<Integer> ptSum : cyclicPart.entrySet()) {
                        int needed = Math.max(spectrumIdx.count(ptSum.getElement()) - nextCandidate.localSpectrum.count(ptSum.getElement()), 0);
                        cyclicScore += ptSum.getCount() - (Math.max(ptSum.getCount() - needed, 0));
                    }
                    if (leaders.isEmpty() || cyclicScore > leaders.get(0).score) {
                        leaders.clear();
                        leaders.add(new Leader(nextCandidate.acids, cyclicScore));
                    } else if (cyclicScore == leaders.get(0).score) {
                        if (nextCandidate.acids.length > leaders.get(0).acids.length) continue;
                        if (nextCandidate.acids.length < leaders.get(0).acids.length) {
                            leaders.clear();
                        }
                        leaders.add(new Leader(nextCandidate.acids, cyclicScore));
                    }
                    // important optimization that removes cheap candidates
                } else if (nextCandidate.mass < spectrum[spectrum.length - 1]) {
                    if (!scoreLen.containsKey(offeredScore)) {
                        scoreLen.put(offeredScore, nextCandidate.acids.length);
                        if (scoreLen.size() > n) {
                            final int outweighed = scoreLen.firstKey();
                            candidates.removeIf(c -> c.score == outweighed);
                            scoreLen.remove(outweighed);
                        }
                    } else scoreLen.put(offeredScore, nextCandidate.acids.length);
                    candidates.add(nextCandidate);
                }
            }
        }
        return leaders.stream().map(x -> x.acids).collect(Collectors.toList());
    }

    public static List<int[]> convolutionCyclopeptideSequencing(final int[] spectrum, final int n, final int m, final int minAcidMass, final int maxAcidMass) {
        Multiset<Integer> convolution = spectrumConvolution(spectrum, minAcidMass, maxAcidMass);
        // https://stackoverflow.com/questions/43806467/java-streams-how-to-do-an-efficient-distinct-and-sort
        Set<Multiset.Entry<Integer>> massFreq = convolution.entrySet();
        int minViableCount = massFreq.stream()
                .map(Multiset.Entry::getCount)
                .sorted(Comparator.reverseOrder())
                .distinct()
                .skip(m - 1)
                .findFirst()
                .orElse(1);
        int[] acidsPool = massFreq.stream()
                .filter(x -> x.getCount() >= minViableCount)
                .mapToInt(Multiset.Entry::getElement)
                .toArray();
        return leaderboardCyclopeptideSequencing(spectrum, acidsPool, n);
    }

    private static Multiset<Integer> spectrumConvolution(final int[] spectrum, final int minAcidMass, final int maxAcidMass) {
        Multiset<Integer> convolution = HashMultiset.create();
        for (int i = 0; i < spectrum.length - 1; i++) {
            for (int j = i + 1; j < spectrum.length; j++) {
                int diff = spectrum[j] - spectrum[i];
                if (diff < minAcidMass) continue;
                if (diff > maxAcidMass) break;
                convolution.add(diff);
            }
        }
        return convolution;
    }

    private static class PeptideCandidate {
        public final int[] acids;
        public final int mass;
        public final Multiset<Integer> localSpectrum;
        public final int score;
        private final int[] lastDiagonal;

        public PeptideCandidate(int[] acids, int mass, Multiset<Integer> localSpectrum, int score, int[] lastDiagonal) {
            this.acids = acids;
            this.localSpectrum = localSpectrum;
            this.mass = mass;
            this.score = score;
            this.lastDiagonal = lastDiagonal;
        }

        public PeptideCandidate() {
            acids = new int[0];
            mass = 0;
            localSpectrum = HashMultiset.create();
            score = 0;
            lastDiagonal = new int[0];
        }

        /**
         * @param newAcid mass of an acid to be added to every spectrum sum
         * @return generated sums
         */
        public int[] deltaSpectrum(int newAcid) {
            int[] delta = new int[lastDiagonal.length + 1];
            for (int i = 0; i < lastDiagonal.length; i++) {
                delta[i] = lastDiagonal[i] + newAcid;
            }
            delta[lastDiagonal.length] = newAcid;
            return delta;
        }

        public PeptideCandidate ascend(int newAcid, int[] recentSubPeptides, int scoreRise) {
            int[] nextAcids = Arrays.copyOf(acids, acids.length + 1);
            nextAcids[acids.length] = newAcid;
            Multiset<Integer> nextLocalSpectrum = HashMultiset.create(localSpectrum);
            for (int subPeptide : recentSubPeptides) nextLocalSpectrum.add(subPeptide);
            return new PeptideCandidate(nextAcids, mass + newAcid, nextLocalSpectrum, score + scoreRise, recentSubPeptides);
        }

        // use Multiset and not array because it allows for utilization of score cache
        public Multiset<Integer> cyclicSpectrumAddition() {
            Multiset<Integer> cyclic = HashMultiset.create();
            int columnFirst = 0;
            for (int i = acids.length - 1; i > 1; i--) {
                int columnAccum = columnFirst += acids[i];
                for (int j = 0; j < i - 1; j++) {
                    cyclic.add(columnAccum += acids[j]);
                }
            }
            return cyclic;
        }

        @Override
        public String toString() {
            return "PeptideCandidate{" +
                    "acids=" + Arrays.toString(acids) +
                    ", localSpectrum=" + localSpectrum +
                    ", score=" + score +
                    ", lastDiagonal=" + Arrays.toString(lastDiagonal) +
                    '}';
        }
    }

    private static class Leader {
        public final int[] acids;
        public final int score;

        public Leader(int[] acids, int score) {
            this.acids = acids;
            this.score = score;
        }
    }
}
