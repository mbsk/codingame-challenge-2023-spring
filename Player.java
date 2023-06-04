import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

/**
 * Auto-generated code below aims at helping you parse
 * the standard input according to the problem statement.
 **/
class Player {

    static final int BRANCH_COUNT = 5;
    static double PI3rd = Math.PI/3;

    static int myAnts=0;
    static int oppAnts=0;
    static Set<Graph> myBases=new HashSet<>();
    static Set<Integer> oppBaseIndex=new HashSet<>();

    static Cell[] cells;
    static int numberOfCells;
    static Set<Integer> resources = new HashSet<>();
    static Set<Integer> eggs = new HashSet<>();
    static Map<Integer, Map<Integer, Integer[]>> cache;
    static Set<Integer> beacons = new HashSet<>();

    static class Cell {
        int type;
        int initialResource;
        int currentResource;
        int[] next = new int[6];
        double x= Double.MIN_VALUE;
        double y= Double.MIN_VALUE;

        int myAnts;
        int oppAnts;
    }

    static class Graph {
        final int rootIdx;
        final Set<Integer> tree = new HashSet<>(); // contains beacon positions
        final List<GraphNode> nodes;
        int segmentsCount = 0;
        long resourcesCount = 0;

		public Graph(int baseIdx, int size) {
            this.rootIdx = baseIdx;
			this.nodes = Stream.generate(() -> new GraphNode()).limit(size).collect(Collectors.toCollection(ArrayList::new));
            IntStream.range(0, size).forEach(i->this.nodes.get(i).idx=i);
            GraphNode rootNode = nodes.get(rootIdx);
            rootNode.distanceToBase = 0;
            rootNode.segment = new ArrayList<>();
            rootNode.segment.add(baseIdx);
            tree.add(baseIdx);
		}

        void updateTree() {
            cleanTree();
            //bestTree(BRANCH_COUNT-segmentsCount);
            //bestTree(BRANCH_COUNT);
            bestTree(0d);
        }

        void cleanTree() {
            tree.clear();
            tree.add(rootIdx);

 
            // cut branches ending with an empty ressource
            /*List<Integer> toBeRemoved = new ArrayList<>();
            tree.stream().filter(idx -> {
                GraphNode node = nodes.get(idx);
                Cell cell = cells[idx];
                //System.err.println("cleanTree1 "+idx + " segment.size="+node.segment.size() + " childs.size="+node.segment.size() + " currentResource="+cell.currentResource);
                // filter leafs (last node in a branch)
                return idx != rootIdx && node.isLeaf() && cell.currentResource == 0;
            }).forEach(idx -> {
                GraphNode leaf = nodes.get(idx);
                int headIdx = leaf.segment.get(0);
                resourcesCount -= leaf.segment.stream().skip(1).filter(i -> cells[i].type>0).count();
                toBeRemoved.addAll(leaf.segment);
                leaf.segment.clear();
                segmentsCount--;
                //System.err.println("cleanTree2 "+idx + " leaf.segment.size="+leaf.segment.size());
                nodes.get(headIdx).childSegments.remove(leaf.segment);
            });
            //System.err.println("toBeRemoved="+toBeRemoved);
            tree.removeAll(toBeRemoved.stream().filter(i->nodes.get(i).segment.size() == 0).collect(Collectors.toList()));
            //System.err.println("tree="+tree);
            updateResourceCount();*/
        }

        /*void bestTree(int depth) {
            Queue<Integer> shortest = null;
            if(depth>0) {
                for(int from: tree) {
                    for (int to : resources) {
                        if(cells[to].currentResource>0 && !tree.contains(to)) { // don't target empty resources and resources already visited by the tree
                            Queue<Integer> path = shortestPathTo(from, to);
                            if(path != null && (shortest == null || path.size() < shortest.size())) {
                                if(Stream.concat(tree.stream(), path.stream()).distinct().count() <= myAnts/(myBases.size()/2.0d)) { // don't bother to go too far 
                                    shortest = path;
                                }
                            }
                        }
                    }
                }
                if(shortest!=null && shortest.size()>1) {
                    //System.err.println("shortest="+shortest);
                    addSegment(new ArrayList<>(shortest));
                    bestTree(depth-1);
                }
            }
        }*/

        // maximize ressource in shorter path (= myAnts)
        void bestTree(double score) {
            System.err.println("bestTree score="+score+" tree.size="+tree.size()+" tree="+tree);
            Deque<Integer> best = null;
            double bestScore = 0d;
            for(int from: tree) {
                for (int to : resources) {
                    if(cells[to].currentResource>0 && !tree.contains(to)) { // don't target empty resources and resources already visited by the tree
                        Deque<Integer> path = shortestPathTo(from, to);
                        //System.err.println("from="+from+" to="+to+" path.size()="+path.size()+" path="+path);
                        if(path != null && tree.size()+path.size()-1 <= myAnts/myBases.size()/(2.0d*myAnts/oppAnts)) {
                            double scoreCandidate = evaluateSegment(path);
                            //System.err.println(" score="+score+" bestScore="+bestScore+" pathScore="+scoreCandidate);
                            if(scoreCandidate > bestScore) {
                                //System.err.println(" =======> BETTER");
                                bestScore = scoreCandidate;
                                best = path;
                            }
                        }
                    }
                }
            }
            if(best!=null && best.size()>1) {
                //System.err.println("shortest="+shortest);
                addSegment(new ArrayList<>(best));
                bestTree(bestScore);
            }
            System.err.println("");
        }

        double evaluateSegment(Collection<Integer> segmentCandidate) { 
            //System.err.println("  evaluateSegment "+segmentCandidate);
            return Double.valueOf((Stream.concat(tree.stream(), segmentCandidate.stream()).filter(i -> cells[i].currentResource>0).count()))/Double.valueOf(segmentCandidate.size());
        }

        void addSegment(List<Integer> segment) {
            //System.err.println("addSegment="+segment);
            tree.addAll(segment);
            int headIdx = segment.get(0);
            GraphNode headNode = nodes.get(headIdx);
            headNode.childSegments.add(segment);
            IntStream.range(1, segment.size()).forEachOrdered(idx -> {
                int nodeIdx = segment.get(idx);
                GraphNode node = nodes.get(nodeIdx);
                node.distanceToBase = headNode.distanceToBase + idx;
                node.segment = segment;
                //System.err.println("  addSegment to "+nodeIdx+" : "+node.distanceToBase);
            });
            segmentsCount++;
            updateResourceCount();
        }

        void updateResourceCount() {
            resourcesCount = tree.stream().filter(i -> cells[i].currentResource > 0).count();
        }

    }

    static class GraphNode {
		List<Integer> segment;
        int distanceToBase = 0;
        int idx;

		final List<List<Integer>> childSegments = new ArrayList<>();

        boolean isLeaf() {
            return segment.size()>0 && idx==segment.get(segment.size()-1) && childSegments.isEmpty();
        }
    }

    static class AstarNode implements Comparable  {

        int idx; // which cell
        int by; // through which cell

        public void computeDistanceTo(int idxEnd) {
            Cell cell1 = cells[idx];
            Cell cell2 = cells[idxEnd];
            this.distanceToEnd = Math.sqrt(Math.pow((cell2.x-cell1.x), 2)+Math.pow((cell2.y-cell1.y), 2));
        }
        double distanceToEnd;
        double weightFromStart;

        double weight() {
            return distanceToEnd+weightFromStart;
        }

        @Override
        public int compareTo(Object o) {
            double diff = weight()-((AstarNode) o).weight();
            return (diff < 0)?-1:(diff>0)?1:0;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;
            AstarNode astarNord = (AstarNode) o;
            return idx == astarNord.idx && by == astarNord.by && Double.compare(astarNord.distanceToEnd, distanceToEnd) == 0 && Double.compare(astarNord.weightFromStart, weightFromStart) == 0;
        }

        @Override
        public int hashCode() {
            return Objects.hash(idx, by, distanceToEnd, weightFromStart);
        }
    }

    static void computePosition(int i) {
        Cell c = cells[i];
        for (int j = 0; j < 6; j++) {
            int inext = c.next[j];
            if (inext >= 0) {
                Cell next = cells[inext];
                if (next.x == Double.MIN_VALUE) {
                    next.x = c.x + Math.cos(PI3rd * j);
                    next.y = c.y + Math.sin(PI3rd * j);
                    computePosition(inext);
                }
            }
        }
    }

    static Deque<Integer> populateCache(int from, int to, AstarNode[] done) {
        int i = to;
        Deque<Integer> result = new ArrayDeque<>();
        // populate the cache
        while(i!=from) {
            result.push(i);
            if(i!=to) {
                cache.get(to).put(i, result.toArray(new Integer[0]));
            }
            i = done[i].by; 
        }
        result.push(from);
        cache.get(to).put(from, result.toArray(new Integer[0]));
        return result;
    }

    static Deque<Integer> shortestPathTo(int from, int to) {
        Integer[] cachedResult = cache.get(to).get(from);
        Deque<Integer> result = new ArrayDeque<>();
        if(cachedResult!=null) {
            Collections.addAll(result, cachedResult);
            //System.err.println("from="+from+" to="+to+" => CACHED "+result);
            return result;
        }
        
        AstarNode[] done = new AstarNode[numberOfCells];
        PriorityQueue<AstarNode> stack = new PriorityQueue<>();
        AstarNode startNode = new AstarNode();
        startNode.idx = from;
        startNode.weightFromStart=0;
        startNode.computeDistanceTo(to);
        stack.add(startNode);
        while (stack.size()>0 && stack.peek().idx != to) {

            // remove first and add it to the "done" list
            AstarNode first = stack.poll();
            done[first.idx] = first;

            // found a ressource, so populate the cache
            if(resources.contains(first.idx)) {
                populateCache(from, first.idx, done);
            }

            // expand first's childs
            Cell cell = cells[first.idx];
            for (int i = 0; i < cell.next.length; i++) {
                int idx = cell.next[i];
                if(idx != -1 && done[idx] == null) {
                    AstarNode node = new AstarNode();
                    node.idx = idx;
                    node.by = first.idx;
                    node.weightFromStart = first.weightFromStart+1;
                    node.computeDistanceTo(to);
                    stack.add(node);
                }
            }
        }
        done[stack.peek().idx] = stack.peek();
        result = populateCache(from, to, done);
        //System.err.println("from="+from+" to="+to+" => "+result);
        return result;
    }

    public static void main(String[] args) {
        Scanner in = new Scanner(System.in);
        numberOfCells = in.nextInt(); // amount of hexagonal cells in this map
        cells = new Cell[numberOfCells];
        for (int i = 0; i < numberOfCells; i++) {
            int type = in.nextInt(); // 0 for empty, 1 for eggs, 2 for crystal
            if(type==1) {
                eggs.add(i);
            }
            if(type==2) {
                resources.add(i);
            }
            Cell cell = new Cell();
            cell.type = type;
            cell.initialResource = in.nextInt(); // the initial amount of eggs/crystals on this cell
            cell.currentResource = cell.initialResource;
            for (int j = 0; j < 6; j++) {
                cell.next[j] = in.nextInt(); // the index of the neighbouring cell for each direction
            }
            cells[i] = cell;
        }
        resources.addAll(eggs);
        cells[0].x = 0;
        cells[0].y = 0;
        computePosition(0);

        int numberOfBases = in.nextInt();
        for (int i = 0; i < numberOfBases; i++) {
            Graph base = new Graph(in.nextInt(), numberOfCells);
            myBases.add(base);
        }
        for (int i = 0; i < numberOfBases; i++) {
            oppBaseIndex.add(in.nextInt());
        }

        cache = new HashMap<>();
        for (int i : resources) {
            cache.put(i, new HashMap<>());
        }

        // precompute shortest distance cache from each base to each ressources
        long before = System.currentTimeMillis();
        for (Graph base : myBases) {
            for (int to : resources) {
                if(cache.get(to).get(base.rootIdx) == null) {
                    shortestPathTo(base.rootIdx, to);
                }
                //System.err.println("from="+from+" to="+to+" dist="+cache.get(to).get(from).length);
            }
        }
        long now = System.currentTimeMillis();
        System.err.println("precompute cache nbCells="+numberOfCells+" nbRes="+resources.size()+" nbBase="+myBases.size()+": "+(now-before)+"ms");


        // game loop
        while (true) {
            myAnts=0;
            oppAnts=0;
            for (int i = 0; i < numberOfCells; i++) {
                cells[i].currentResource = in.nextInt(); // the current amount of eggs/crystals on this cell
                cells[i].myAnts = in.nextInt(); // the amount of your ants on this cell
                myAnts+=cells[i].myAnts;
                cells[i].oppAnts = in.nextInt(); // the amount of opponent ants on this cell
                oppAnts+=cells[i].oppAnts;
            }

            beacons.clear();
            for (Graph graph : myBases) {
                graph.updateTree();
                beacons.addAll(graph.tree);
            }
            String cmd="";
            for (Integer cellIdx : beacons) {
                cmd+="BEACON "+cellIdx+" 1;";
            }
            System.out.println(cmd);



            // Write an action using System.out.println()
            // To debug: System.err.println("Debug messages...");


            // WAIT | LINE <sourceIdx> <targetIdx> <strength> | BEACON <cellIdx> <strength> | MESSAGE <text>
            //System.out.println("WAIT");
        }
    }
}
