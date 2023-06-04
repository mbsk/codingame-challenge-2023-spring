import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

class Player {

    static long chronoTop = 0L;

    final static double PI3rd = Math.PI/3;

    static int numberOfCells;       // board's cells count
    static Cell[] cells;            // board of the game, an array of cells indexed by cell's idx
    static Set<Integer> resources = new HashSet<>();    // all resource (crystal and eggs) of the game indexed by cell's idx

    static Opponent opp = new Opponent();       // holds opponent data
    static Strategy strategy = new Strategy();  // holds chain building strategy

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

    /** Holds opponent's data */
    static class Opponent {
        int ants=0;        // counter of opponent ants
        int score=0;
        final Set<Integer> bases = new HashSet<>();     // list of opponent bases indexed by cell's idx
        final Set<Integer> occupied = new HashSet<>();    // list of opponent occupied cells indexed by cell's idx

        void resetCounters() {
            ants = 0;
            occupied.clear();
        }
        public void addAnts(int i, int oppAnts) {
            if(oppAnts > 0){
                occupied.add(i);
                ants+=oppAnts;
            }
        }
    }

    /** Represents an node for the Astar shotest path finding algorithm */
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

        public int compareTo(Object o) {
            double diff = weight()-((AstarNode) o).weight();
            return (diff < 0)?-1:(diff>0)?1:0;
        }

        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;
            AstarNode astarNord = (AstarNode) o;
            return idx == astarNord.idx && by == astarNord.by && Double.compare(astarNord.distanceToEnd, distanceToEnd) == 0 && Double.compare(astarNord.weightFromStart, weightFromStart) == 0;
        }

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

    /** Holds the strategy to create a chain of beacons to harvest resources */
    static class Strategy {
        int initialEggCount = 0;
        int initialCrystalCount = 0;

        final Set<Integer> bases = new HashSet<>(); // our bases cell's idx
        final Set<Integer> chain = new HashSet<>();  // an array of beacon positions (cell's idx)
        final Set<Integer> harvestingResources = new HashSet<>(); // resources currently harvested by the chain

        int ants=0;               // counter of own ants
        int score=0;
        final Set<Integer> occupied = new HashSet<>();    // list of opponent occupied cells indexed by cell's idx
        int eggsLeft=0;             // counter of remaining eggsCount
        int crystalsLeft =0;        // counter of remaining resources

        public void addAnts(int i, int oppAnts) {
            if(oppAnts > 0){
                occupied.add(i);
                ants+=oppAnts;
            }
        }

        void resetCounters() {
            occupied.clear();
            ants = 0;
            eggsLeft = 0;
            crystalsLeft = 0;
        }

        // called on each turn to update the tree
        void updateChain() {
            System.err.println(System.currentTimeMillis()-chronoTop+" updateChain()");
            // clear the chain only if one of the harvested resources became empty
            if(harvestingResources.size() ==0 || chain.stream().filter(i-> cells[i].type>0 && cells[i].currentResource>0).count() != harvestingResources.size()) {
                System.err.println(System.currentTimeMillis()-chronoTop+" updateChain() : need refresh");
                chain.clear();
                bases.stream().forEach(i -> chain.add(i));
                addBranch(0d, 10);
                System.err.println(System.currentTimeMillis()-chronoTop+" finish add branch");
                harvestingResources.clear();
                chain.stream().filter(i-> cells[i].type>0 && cells[i].currentResource>0).forEach(i -> harvestingResources.add(i));
            }
        }

        /*  find a new branch from all current tree position which maximize resources in shorter path */
        void addBranch(double score, int iteration) {
            if(iteration==0) return;
            //System.err.println("bestTree score="+score+" tree.size="+tree.size()+" tree="+tree);
            Deque<Integer> newBranch = null;
            double newBranchScore = 0d;
            for(int from: chain) {
                for (int to: resources) {
                    // loop : calculate shorter path from each current tree position to all resources
                    //System.err.println("res : "+to+" "+cells[to].currentResource);
                    if(cells[to].currentResource>0 && !chain.contains(to)) { // don't target empty resources and resources already visited by the tree
                        double maxTreeSize = ants;
                        Deque<Integer> pathCandidate = shortestPathTo(from, to);
                        //System.err.println("from="+from+" to="+to+" path.size()="+path.size()+" path="+path);
                        if(pathCandidate != null && Double.valueOf(chain.size()+pathCandidate.size()-1) <= maxTreeSize) {
                            double scoreCandidate = evaluateSegment(pathCandidate);
                            //System.err.println(" score="+score+" newBranchScore="+newBranchScore+" pathScore="+scoreCandidate);
                            if(scoreCandidate > newBranchScore) {
                                //System.err.println(" =======> BETTER");
                                newBranchScore = scoreCandidate;
                                newBranch = pathCandidate;
                            }
                        }
                    }
                }
            }
            // if a new branch was found, add it to the graph, and find a new one
            if(newBranch!=null && newBranch.size()>1) {
                //System.err.println("shortest="+shortest);
                System.err.println(System.currentTimeMillis()-chronoTop+" add branch");
                addSegment(new ArrayList<>(newBranch));
                addBranch(newBranchScore, iteration-1);
            }
        }

        // a segment score
        double evaluateSegment(Collection<Integer> segmentCandidate) {
            List<Integer> resourcesLeft = resources.stream().filter(i -> cells[i].type==2 && cells[i].currentResource>0).collect(Collectors.toList());
            if(resourcesLeft.size()==1 && !isReachable(resourcesLeft.get(0)) && !segmentCandidate.contains(resourcesLeft.get(0))) return 0; // if there is only 1 reachable resource, we should rush for it
            final double[] score = { 0d, 0d}; // an "effective static" wrapper to be used in lambda score[0]:"sum of cell's score" score[1]:"sum of cells"
            double count = Double.valueOf(Stream.concat(chain.stream(), segmentCandidate.stream()).distinct().count());
            score[1] = Math.pow(count,3);
            Stream.concat(chain.stream(), segmentCandidate.stream()).filter(i -> cells[i].currentResource>0).distinct().forEach(i->{
                if(cells[i].type == 1) {
                    // it's an EGG
                    score[0]++;
                } else if (cells[i].type ==2) {
                    // it's a crystal
                    score[0]++;
                    score[0]+= cells[i].currentResource/ crystalsLeft; // prioritize most valueable resource
                } else {
                    score[1]*=2; // penalize empty cells
                }
            });
            return score[0]/score[1];
        }

        boolean isReachable(int to) {
            boolean result = false;
            for (int from:bases){
                result = result || shortestPathTo(from, to).size() <= ants;
            }
            return result;
        }

        void addSegment(List<Integer> segment) {
            //System.err.println("addSegment="+segment);
            chain.addAll(segment);
        }

        public void addInitialResources(int type, int count) {
            if(type==1) {
                initialEggCount+=count;
            } else if (type==2) {
                initialCrystalCount+=count;
            }
        }
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
            strategy.addInitialResources(cell.type, cell.initialResource);
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
            strategy.bases.add(in.nextInt());
        }
        for (int i = 0; i < numberOfBases; i++) {
            opp.bases.add(in.nextInt());
        }

        shortestPathCache = new HashMap<>();
        for (int i : resources) {
            shortestPathCache.put(i, new HashMap<>());
        }

        // precompute shortest distance cache from each base to each resources
        chronoTop = System.currentTimeMillis();
        for (int from : strategy.bases) {
            for (int to : resources) {
                if(shortestPathCache.get(to).get(from) == null) {
                    shortestPathTo(from, to);
                }
            }
        }
        long now = System.currentTimeMillis();
        System.err.println("precompute cache nbCells="+numberOfCells+" nbRes="+resources.size()+" nbBase="+ strategy.bases.size()+": "+(now-chronoTop)+"ms");

        Set<Integer> beacons = new HashSet<>(); // holds a list of cells where to put beacons (recalculated every turn)
        // game loop
        while (true) {
            strategy.score = in.nextInt();
            opp.score = in.nextInt();
            
            strategy.resetCounters();
            opp.resetCounters();
            System.err.println(System.currentTimeMillis()-chronoTop+" RESET COUNTERS DONE");
            for (int i = 0; i < numberOfCells; i++) {
                cells[i].currentResource = in.nextInt(); // the current amount of eggs/crystals on this cell
                cells[i].myAnts = in.nextInt(); // the amount of your ants on this cell
                strategy.addAnts(i, cells[i].myAnts);
                cells[i].oppAnts = in.nextInt(); // the amount of opponent ants on this cell
                opp.addAnts(i, cells[i].oppAnts);

                if(cells[i].type==1) {
                    strategy.eggsLeft +=cells[i].currentResource;
                } else if(cells[i].type==2) {
                    strategy.crystalsLeft +=cells[i].currentResource;
                }
            }
            chronoTop = System.currentTimeMillis();
            System.err.println(System.currentTimeMillis()-chronoTop+" START");

            //System.err.println("my.ants="+strategy.ants+" my.occupied="+strategy.occupied.size()+" opp.ant="+opp.ants+" opp.occupied="+opp.occupied.size()+" eggs="+ strategy.eggsLeft +" crystals="+ strategy.crystalsLeft);

            // WAIT | LINE <sourceIdx> <targetIdx> <strength> | BEACON <cellIdx> <strength> | MESSAGE <text>
            beacons.clear();
            strategy.updateChain();
            beacons.addAll(strategy.chain);
            String cmd="";
            for (Integer cellIdx : beacons) {
                cmd+="BEACON "+cellIdx+" 1;";
            }
            System.err.println(System.currentTimeMillis()-chronoTop+" WRITING BEACON");
            System.out.println(cmd);
            System.err.println(System.currentTimeMillis()-chronoTop+" FINISHED loop");

            // Write an action using System.out.println()
            // To debug: System.err.println("Debug messages...");

        }
    }
}
