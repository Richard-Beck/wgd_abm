package metastasis;

import HAL.GridsAndAgents.AgentGrid2D;
import HAL.GridsAndAgents.AgentSQ2D;
import HAL.Gui.OpenGL2DWindow;
import HAL.Tools.FileIO;
import HAL.Rand;
import HAL.Gui.PlotLine;
import HAL.Gui.PlotWindow;
import HAL.Util;

import java.awt.event.KeyListener;
import java.io.IOException;
import java.util.*;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static org.lwjgl.glfw.GLFW.*;



//cells grow and mutate
public class metastasis implements Viz.KeyListener {
    Rand rn=new Rand(1);
    FileIO outputFile=null;
    CTCs C=new CTCs();

    Viz vis;
    int maxGrids=3;
    int x=200;
    int y=200;

    ArrayList<tissueSiteGrid> sites = new ArrayList<>();
    int window=0;
    public metastasis() {

    }

    @Override
    public void onKeyPress(int key) {
        window+=1;
        window=window%sites.size();
        System.out.println("switching to window " + (window+1) + " of " + sites.size());
    }

    public void addSite(CellEx c){
        tissueSiteGrid met = new tissueSiteGrid(x,y,vis,rn);
        met.InitTumor(-1);
        met.M = new gridManager(x,y);
        int coord = met.RandomAgent(rn).Isq();
        CellEx d = met.NewAgentSQ(coord);
        d.Init(c.type,c.genome);
        met.tumorPop=1;
        sites.add(met);
    }
    public static void main(String[]args) throws IOException {

        IO io = new IO("metastasis/trial_output","metastasis/parameters.csv");

        int x=200,y=200;
        // NOTE had to modify the HAL library to make context current within update to allow 2 windows
        metastasis sim = new metastasis();
        sim.vis = new Viz("Primary", 9*x, 3*y, 3*x, y);
        sim.vis.addKeyListener(sim);
        tissueSiteGrid primary = new tissueSiteGrid(x,y,sim.vis,sim.rn);
        primary.InitTumor(5);
        io.initializeCellTypes();
        primary.M = new gridManager(x,y);
        sim.sites.add(primary);

        for(int tick = 0; tick<20000; tick++){
            if(sim.vis.IsClosed()) break;
            ArrayList<CellEx> ex= sim.C.stepCTCs(sim.rn);
            int n2add = min(max(sim.maxGrids-sim.sites.size(),0),ex.size());

            while(n2add>0){
                sim.addSite(ex.get(n2add-1));
                n2add--;
            }

            for(int i = 0; i<sim.sites.size(); i++){
                sim.sites.get(i).stepTumor(tick,sim.C,sim.window==i);
            }

            Iterator<tissueSiteGrid> iterator = sim.sites.iterator();
            while (iterator.hasNext()) {
                int value = iterator.next().tumorPop;
                if (value == 0) { // Remove failed mets.
                    iterator.remove();
                }
            }

        }

        io.writeAllCells(sim.sites);
    }
}