import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

class Player {

    final static double PI3rd = Math.PI/3;

    static int myAnts=0;    // counter of own ants
    static int oppAnts=0;   // counter of opponent ants
    static int eggsCount=0; // counter of initial eggsCount

    static Set<Graph> myBases=new HashSet<>();          // list of own bases
    static Set<Integer> oppBaseIndex=new HashSet<>();   // list of opponent bases

    static Cell[] cells;        // board of the game, an array of cells indexed by cell's idx
    static int numberOfCells;   // board's cells count
    static Set<Integer> resources = new HashSet<>();    // all resource (crystal and eggs) of the game indexed by cell's idx

    static Map<Integer, Map<Integer, Integer[]>> shortestPathCache; // holds the cache of shorter path to each resource (to, from)

    /** Represents a cell of the board */
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

    /** Represents a tree of beacons from a base */
    static class Graph {
        final int rootIdx; // cell's idx of the base
        final Set<Integer> tree = new HashSet<>(); // an array of beacon positions (cell's idx)
        final List<GraphNode> nodes; // graph of beacon's tree (indexed by cell's idx)
        long resourcesCount = 0;    // count of resources covered by the tree

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

        // called on each turn to update the tree
        void updateTree() {
            cleanTree();
            //bestTree(BRANCH_COUNT-segmentsCount);
            //bestTree(BRANCH_COUNT);
            bestTree(0d);
        }

        void cleanTree() {
            // STRATEGY 1 : clean the tree on each turn
            tree.clear();
            tree.add(rootIdx);

            // STRATEGY 2 : only cut "dead" branches (= the one who leads to an empty ressource leaf) : seems to be less efficient
            // cut branches ending with an empty resource
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

        /*  find a new branch from all current tree position which maximize resources in shorter path */
        void bestTree(double score) {
            //System.err.println("bestTree score="+score+" tree.size="+tree.size()+" tree="+tree);
            Deque<Integer> newBranch = null;
            double newBranchScore = 0d;
            for(int from: tree) {
                for (int to : resources) {
                    // loop : calculate shorter path from each current tree position to all resources
                    if(cells[to].currentResource>0 && !tree.contains(to)) { // don't target empty resources and resources already visited by the tree

                        Deque<Integer> path = shortestPathTo(from, to);
                        //System.err.println("from="+from+" to="+to+" path.size()="+path.size()+" path="+path);
                        if(path != null && tree.size()+path.size()-1 <= myAnts/myBases.size()/(2.0d*myAnts/oppAnts)) {
                            double scoreCandidate = evaluateSegment(path);
                            //System.err.println(" score="+score+" newBranchScore="+newBranchScore+" pathScore="+scoreCandidate);
                            if(scoreCandidate > newBranchScore) {
                                //System.err.println(" =======> BETTER");
                                newBranchScore = scoreCandidate;
                                newBranch = path;
                            }
                        }
                    }
                }
            }
            // if a new branch was found, add it to the graph, and find a new one
            if(newBranch!=null && newBranch.size()>1) {
                //System.err.println("shortest="+shortest);
                addSegment(new ArrayList<>(newBranch));
                bestTree(newBranchScore);
            }
            //System.err.println("");
        }

        // a segment score
        double evaluateSegment(Collection<Integer> segmentCandidate) {
            return Double.valueOf(Stream.concat(tree.stream(), segmentCandidate.stream()).filter(i -> cells[i].currentResource>0).distinct().count())/Double.valueOf(Stream.concat(tree.stream(), segmentCandidate.stream()).distinct().count());
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
            updateResourceCount();
        }

        void updateResourceCount() {
            resourcesCount = tree.stream().filter(i -> cells[i].currentResource > 0).count();
        }

    }

    /** Represents a node in the tree graph */
    static class GraphNode {
        int idx;                // cell's idx
        List<Integer> segment;  // segment of the tree to which this node belongs (an array of cell's idx)
        int distanceToBase = 0; // node's count to base
        final List<List<Integer>> childSegments = new ArrayList<>();    // array of child segments if this node is a bifurcation in the tree

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

    /** calcaulate Cartesian's (x,y) coordinates of each cell of the board */
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

    /** Populate cache whith all intermediate path from 'from' to 'to' */
    static Deque<Integer> populateCache(int from, int to, AstarNode[] path) {
        int i = to;
        Deque<Integer> result = new ArrayDeque<>();
        while(i!=from) {
            result.push(i);
            if(i!=to) {
                shortestPathCache.get(to).put(i, result.toArray(new Integer[0]));
            }
            i = path[i].by;
        }
        result.push(from);
        shortestPathCache.get(to).put(from, result.toArray(new Integer[0]));
        return result;
    }

    /** Astar implementation of shortestPath */
    static Deque<Integer> shortestPathTo(int from, int to) {
        Integer[] cachedResult = shortestPathCache.get(to).get(from);
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

            // found a resource, so populate the cache
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
        numberOfCells = in.nextInt();

        // build
        cells = new Cell[numberOfCells];
        for (int i = 0; i < numberOfCells; i++) {
            int type = in.nextInt(); // 0 for empty, 1 for eggs, 2 for crystal
            if(type>0) {
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

        shortestPathCache = new HashMap<>();
        for (int i : resources) {
            shortestPathCache.put(i, new HashMap<>());
        }

        // precompute shortest distance cache from each base to each resources
        long before = System.currentTimeMillis();
        for (Graph base : myBases) {
            for (int to : resources) {
                if(shortestPathCache.get(to).get(base.rootIdx) == null) {
                    shortestPathTo(base.rootIdx, to);
                }
            }
        }
        long now = System.currentTimeMillis();
        System.err.println("precompute cache nbCells="+numberOfCells+" nbRes="+resources.size()+" nbBase="+myBases.size()+": "+(now-before)+"ms");

        Set<Integer> beacons = new HashSet<>(); // holds a list of cells where to put beacons (recalculated every turn)
        // game loop
        while (true) {
            for (int i = 0; i < numberOfCells; i++) {
                cells[i].currentResource = in.nextInt(); // the current amount of eggs/crystals on this cell
                cells[i].myAnts = in.nextInt(); // the amount of your ants on this cell
                myAnts+=cells[i].myAnts;
                cells[i].oppAnts = in.nextInt(); // the amount of opponent ants on this cell
                oppAnts+=cells[i].oppAnts;
            }

            // WAIT | LINE <sourceIdx> <targetIdx> <strength> | BEACON <cellIdx> <strength> | MESSAGE <text>
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

        }
    }
}

