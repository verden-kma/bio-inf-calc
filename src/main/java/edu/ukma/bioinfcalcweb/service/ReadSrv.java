package edu.ukma.bioinfcalcweb.service;

import com.google.common.collect.Lists;
import org.apache.commons.collections4.ListUtils;
import org.springframework.data.util.Pair;
import org.springframework.util.SerializationUtils;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class ReadSrv {
    public static String patternComposition(String[] patterns) {
        Map<String, List<String>> dbGraph = deBruign(patterns);
        Pair<Map<Integer, String>, Map<Integer, Deque<Integer>>> transformed = transformGraph(dbGraph);
        Map<Integer, String> dec = transformed.getFirst();
        Map<Integer, Deque<Integer>> intGraph = transformed.getSecond();
        List<Integer> path = eulerianPath(intGraph);
        List<String> patternPath = path.stream()
                .map(dec::get)
                .collect(Collectors.toList());
        StringBuilder resDna = new StringBuilder();
        patternPath.forEach(p -> resDna.append(p.charAt(0)));
        resDna.append(patternPath.get(patternPath.size() - 1).substring(1));
        return resDna.toString();
    }

    public static List<String> findContigs(String[] dna) {
        Map<String, List<String>> graph = deBruign(dna);
        Pair<Map<Integer, String>, Map<Integer, Deque<Integer>>> transformed = transformGraph(graph);
        Map<Integer, String> decoder = transformed.getFirst();
        Map<Integer, Deque<Integer>> numGraph = transformed.getSecond();
        Set<Integer> seeds = branchingNodes(numGraph);
        List<List<Integer>> numContigs = new ArrayList<>();
        for (Integer seed : seeds) {
            Deque<Integer> seedAdj = numGraph.get(seed);
            while (!seedAdj.isEmpty()) {
                List<Integer> contig = new ArrayList<>(Arrays.asList(seed));
                Integer v = seedAdj.pop();
                while (!seeds.contains(v) && numGraph.containsKey(v) && numGraph.get(v).size() > 0) {
                    contig.add(v);
                    v = numGraph.get(v).pop();
                }
                contig.add(v);
                numContigs.add(contig);
            }
        }
        List<String> contigs = new ArrayList<>();
        for (List<Integer> numContig : numContigs) {
            List<String> splitContig = numContig.stream().map(decoder::get).collect(Collectors.toList());
            List<String> contigBatch = new ArrayList<>(Arrays.asList(getPrefix(splitContig.get(0))));
            contigBatch.addAll(splitContig.stream().map(x -> String.valueOf(x.charAt(x.length() - 1))).collect(Collectors.toList()));
            contigs.add(String.join("", contigBatch));
        }
        return contigs;
    }

    private static Map<String, List<String>> deBruign(String[] patterns) {
        Map<String, List<String>> prefMap = new HashMap<>();
        for (String pattern : patterns) {
            String prefix = getPrefix(pattern);
            prefMap.putIfAbsent(prefix, new ArrayList<>());
            prefMap.get(prefix).add(pattern);
        }
        Map<String, List<String>> adjList = new HashMap<>();
        for (Map.Entry<String, List<String>> prefPat : prefMap.entrySet()) {
            adjList.put(prefPat.getKey(), new ArrayList<>(prefPat.getValue().size()));
            for (String pattern : prefPat.getValue()) {
                adjList.get(prefPat.getKey()).add(getSuffix(pattern));
            }
        }
        return adjList;
    }

    private static Map<Pair<String, String>, List<Pair<String, String>>> pairedDeBruign(List<Pair<String, String>> patterns) {
        Map<Pair<String, String>, List<Pair<String, String>>> prefMap = new HashMap<>();
        for (Pair<String, String> readPair : patterns) {
            Pair<String, String> nextPrefix = Pair.of(getPrefix(readPair.getFirst()), getPrefix(readPair.getSecond()));
            prefMap.putIfAbsent(nextPrefix, new ArrayList<>());
            prefMap.get(nextPrefix).add(readPair);
        }
        Map<Pair<String, String>, List<Pair<String, String>>> adjList = new HashMap<>();
        for (Pair<String, String> prefPair : prefMap.keySet()) {
            adjList.put(prefPair, new ArrayList<>());
            for (Pair<String, String> patternPair : prefMap.get(prefPair))
                adjList.get(prefPair).add(Pair.of(getSuffix(patternPair.getFirst()), getSuffix(patternPair.getSecond())));
        }
        return adjList;
    }

    private static String getPrefix(String str) {
        return str.substring(0, str.length() - 1);
    }

    private static String getSuffix(String str) {
        return str.substring(1);
    }

    private static Pair<Map<Integer, String>, Map<Integer, Deque<Integer>>> transformGraph(Map<String, List<String>> graph) {
        int indexAccum = 0;
        Map<String, Integer> mapper = new HashMap<>();
        Map<Integer, String> remapper = new HashMap<>();
        Map<Integer, Deque<Integer>> outGraph = new HashMap<>();
        for (Map.Entry<String, List<String>> node : graph.entrySet()) {
            if (!mapper.containsKey(node.getKey())) {
                mapper.put(node.getKey(), ++indexAccum);
                remapper.put(indexAccum, node.getKey());
            }
            outGraph.put(mapper.get(node.getKey()), new ArrayDeque<>());
            for (String adjV : node.getValue()) {
                if (!mapper.containsKey(adjV)) {
                    mapper.put(adjV, ++indexAccum);
                    remapper.put(indexAccum, adjV);
                }
                outGraph.get(mapper.get(node.getKey())).add(mapper.get(adjV));
            }
        }
        return Pair.of(remapper, outGraph);
    }

    private static <T> Pair<Map<Integer, T>, Map<Integer, Deque<Integer>>> transformGraph2(Map<T, List<T>> graph) {
        int indexAccum = 0;
        Map<T, Integer> mapper = new HashMap<>();
        Map<Integer, T> remapper = new HashMap<>();
        Map<Integer, Deque<Integer>> outGraph = new HashMap<>();
        for (Map.Entry<T, List<T>> node : graph.entrySet()) {
            if (!mapper.containsKey(node.getKey())) {
                mapper.put(node.getKey(), ++indexAccum);
                remapper.put(indexAccum, node.getKey());
            }
            outGraph.put(mapper.get(node.getKey()), new ArrayDeque<>());
            for (T adjV : node.getValue()) {
                if (!mapper.containsKey(adjV)) {
                    mapper.put(adjV, ++indexAccum);
                    remapper.put(indexAccum, adjV);
                }
                outGraph.get(mapper.get(node.getKey())).add(mapper.get(adjV));
            }
        }
        return Pair.of(remapper, outGraph);
    }

    private static List<Integer> eulerianPath(Map<Integer, Deque<Integer>> intGraph) {
        Optional<Integer> maybeMax = Stream.concat(intGraph.keySet().stream(),
                intGraph.values().stream().flatMap(Collection::stream)).max(Integer::compareTo);
        if (!maybeMax.isPresent()) return Collections.emptyList();
        int[] dDegrees = new int[maybeMax.get() + 1];
        for (Map.Entry<Integer, Deque<Integer>> node : intGraph.entrySet()) {
            dDegrees[node.getKey()] += node.getValue().size();
            for (Integer v : node.getValue())
                dDegrees[v]--;
        }
        int starter = -1;
        for (int i = 0; i < dDegrees.length; i++) {
            if (dDegrees[i] == 1) starter = i;
        }
        if (starter == -1) return Collections.emptyList();
        Deque<Map.Entry<Integer, Deque<Integer>>> stack = new ArrayDeque<>();
        stack.add(new AbstractMap.SimpleEntry<>(starter, intGraph.get(starter)));
        List<Integer> cycle = new ArrayList<>();
        while (!stack.isEmpty()) {
            Map.Entry<Integer, Deque<Integer>> node = stack.removeLast();
            Integer v = node.getKey();
            Deque<Integer> adj = node.getValue();
            while (!adj.isEmpty()) {
                stack.add(new AbstractMap.SimpleEntry<>(v, adj));
                v = adj.removeFirst();
                adj = intGraph.getOrDefault(v, new ArrayDeque<>());
            }
            cycle.add(v);
        }
        return Lists.reverse(cycle);
    }

    private static List<List<Integer>> allEulerianPaths(Map<Integer, Deque<Integer>> graph, Deque<Integer> stack) {
        List<Integer> singlePath = new ArrayList<>();
        List<List<Integer>> paths = new ArrayList<>();
        while (!stack.isEmpty()) {
            Integer v = stack.removeLast();
            Deque<Integer> adj = graph.getOrDefault(v, new ArrayDeque<>());
            while (!adj.isEmpty()) {
                stack.add(v);
                if (adj.size() == 1) {
                    v = adj.removeFirst();
                    adj = graph.getOrDefault(v, new ArrayDeque<>());
                } else {
                    for (Integer adjV : adj) {
                        @SuppressWarnings("unchecked")
                        Map<Integer, Deque<Integer>> localGraph = Objects.requireNonNull(graph.getClass()
                                .cast(SerializationUtils.deserialize(SerializationUtils.serialize(graph))));
                        localGraph.get(v).remove(adjV);
                        @SuppressWarnings("unchecked")
                        Deque<Integer> localStack = Objects.requireNonNull(stack.getClass()
                                .cast(SerializationUtils.deserialize(SerializationUtils.serialize(stack))));
                        localStack.add(adjV);
                        for (List<Integer> pathTail : allEulerianPaths(localGraph, localStack))
                            paths.add(ListUtils.union(singlePath, pathTail));
                        return paths;
                    }
                }
            }
            singlePath.add(v);
        }
        return Collections.singletonList(singlePath);
    }

    private static Set<Integer> branchingNodes(Map<Integer, Deque<Integer>> numGraph) {
        Map<Integer, List<Integer>> adjIn = new HashMap<>();
        for (Map.Entry<Integer, Deque<Integer>> adjEntry : numGraph.entrySet())
            for (Integer adjV : adjEntry.getValue()) {
                if (!adjIn.containsKey(adjV))
                    adjIn.put(adjV, new ArrayList<>());
                adjIn.get(adjV).add(adjEntry.getKey());
            }
        Set<Integer> res = new HashSet<>();
        for (Integer v : numGraph.keySet())
            if (numGraph.get(v).size() > 1) res.add(v);
            else if (numGraph.get(v).size() == 1 && !adjIn.containsKey(v) || adjIn.get(v).size() > 1) res.add(v);
        return res;
    }

    public static String pairedPatternReconstruction(List<Pair<String, String>> vertices, int d) {
        int k = vertices.get(0).getFirst().length() - 1;
        Map<Pair<String, String>, List<Pair<String, String>>> originalGraph = pairedDeBruign(vertices);
        Pair<Map<Integer, Pair<String, String>>, Map<Integer, Deque<Integer>>> transformation = transformGraph2(originalGraph);
        Map<Integer, Pair<String, String>> remapper = transformation.getFirst();
        Map<Integer, Deque<Integer>> mappedGraph = transformation.getSecond();
        List<List<Integer>> paths = allEulerianPaths(mappedGraph, new ArrayDeque<>(Collections.singletonList(peekStart(mappedGraph))))
                .stream().distinct().collect(Collectors.toList());
        paths.forEach(Collections::reverse);
        List<List<Pair<String, String>>> decodedPaths = new ArrayList<>(paths.size());
        for (List<Integer> encPath : paths)
            decodedPaths.add(encPath.stream().map(remapper::get).collect(Collectors.toList()));
        List<Pair<String, String>> validPath = null;
        for (List<Pair<String, String>> path : decodedPaths) {
            boolean pathFailed = false;
            for (int i = 0; i < path.size(); i++) {
                String fstRead = path.get(i).getFirst();
                String sndRead = path.get(i).getSecond();
                if (i == 0) {
                    for (int j = 0; j < k - 1; j++) {
                        pathFailed = !checkColumn(fstRead, sndRead, i, j, k, d, path);
                        if (pathFailed) break;
                    }
                    if (pathFailed) break;
                }
                pathFailed = !checkColumn(fstRead, sndRead, i, k - 1, k, d, path);
                if (pathFailed) break;
            }
            if (pathFailed) continue;
            validPath = path;
        }
        final StringBuilder genomePrefix = new StringBuilder();
        if (validPath == null) return null;
        validPath.forEach(node -> genomePrefix.append(node.getFirst().charAt(0)));
        genomePrefix.append(validPath.get(validPath.size() - 1).getFirst().substring(1));
        for (int i = d + 1; i > 0; i--)
            genomePrefix.append(validPath.get(validPath.size() - 1 - i).getSecond().charAt(0));
        genomePrefix.append(validPath.get(validPath.size() - 1).getSecond());
        return genomePrefix.toString();
    }

    private static boolean checkColumn(String fstRead, String sndRead, int i, int j, int k, int d, List<Pair<String, String>> path) {
        for (int p = 1; p < k; p++) {
            // current index is > than last index of a read AND read pair that should have contained curr acid is not out of scope
            // AND 'right angle' assertion for 1st read failed
            if (j >= p && i + p < path.size() && fstRead.charAt(j) != path.get(i + p).getFirst().charAt(j - p))
                return false;
            int absPos = k + d + j - p;
            Boolean pair2Proj = null; // is in 'd' split
            if (absPos < k) // is in 1st read
                pair2Proj = true;
            else if (absPos >= k + d) // is in 2nd read
                pair2Proj = false;
            // proj in not in split AND right angle assertion is possible AND right angle assertion for 2nd read failed
            if (pair2Proj != null && i + p < path.size()
                    && sndRead.charAt(j) != (pair2Proj ? path.get(i + p).getFirst().charAt(absPos)
                    : path.get(i + p).getSecond().charAt(absPos - k - d)))
                return false;
        }
        return true;
    }

    private static Integer peekStart(Map<Integer, Deque<Integer>> mappedGraph) {
        Optional<Integer> maybeMax = Stream.concat(mappedGraph.keySet().stream(),
                mappedGraph.values().stream().flatMap(Collection::stream)).max(Integer::compareTo);
        if (!maybeMax.isPresent()) throw new NoSuchElementException("Could not find maximum vertex label.");
        int[] dDegrees = new int[maybeMax.get() + 1];
        for (Integer v : mappedGraph.keySet()) {
            Deque<Integer> adj = mappedGraph.get(v);
            dDegrees[v] += adj.size();
            for (Integer vx : adj)
                dDegrees[vx]--;
        }
        Integer starter = null;
        for (int i = 0; i < dDegrees.length; i++)
            if (dDegrees[i] == 1) starter = i;
        if (starter == null) throw new NoSuchElementException("Could not find a peek.");
        return starter;
    }
}
