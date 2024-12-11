package metastasis;

import HAL.Gui.OpenGL2DWindow;
import HAL.Tools.FileIO;
import HAL.Util;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;

import static HAL.Util.RGB256;

public class IO {

    int imageWriteFrequency=100;
    String outputPath;
    Boolean disable_output=false;
    File outputFolder;
    File imageFolder;
    FileIO summary=null;
    IO(String path2output,String path2Params) throws IOException {
        File outputHome = new File(path2output);
        File[] folderList = outputHome.listFiles();
        outputPath = outputHome.getAbsolutePath();// + "/run_" + folderList.length;
        outputFolder = new File(outputPath);
        outputFolder.mkdirs();
        //make a copy of the parameter file in the output directory:
        Path paramSource = Paths.get(path2Params);
        Path paramDest = Paths.get(outputPath + "/params.txt");
        Files.copy(paramSource, paramDest, StandardCopyOption.REPLACE_EXISTING);

        // create subfolders to store any images.
        imageFolder = new File(outputPath+"/images");
        imageFolder.mkdirs();

        summary = new FileIO(outputPath + "/summary.csv", "w");
        summary.Write("time,site,tumor,wgd,vessel\n");
    }

    IO(){
        disable_output = true;
    };

    public void initializeCellTypes()  {
        FileIO reader = new FileIO(outputPath + "/params.txt","r");
        String line = reader.ReadLine(); // Skip header

        while ((line = reader.ReadLine()) != null) {
            String[] values = line.split(",");
            String type = values[0];
            double MMPGen = Double.parseDouble(values[1]);
            double ECMDeg = Double.parseDouble(values[2]);
            double HaptoCoef= Double.parseDouble(values[3]);
            double pMove= Double.parseDouble(values[4]);
            int r = Integer.parseInt(values[5]);
            int g = Integer.parseInt(values[6]);
            int b = Integer.parseInt(values[7]);
            double space = Double.parseDouble(values[8]);
            int color = Util.RGB(r, g, b);
            // Initialize the corresponding enum value
            CellEx.CellType.valueOf(type).setParams(MMPGen, ECMDeg, HaptoCoef,pMove,color, space);

        }

        reader.Close();
    }

    public void Draw(int time, OpenGL2DWindow vis){
        if(disable_output) return;
        if(time%imageWriteFrequency!=0) return;
        if(vis==null) return;
        vis.ToJPG(imageFolder.getAbsolutePath()+"/"+time+".jpg");
    }

    public void writeAllCells(ArrayList<tissueSiteGrid> dataWriters) throws IOException {
        FileIO writer = new FileIO(imageFolder + "/test.csv","w");
        int gridNo = 0;
        for(tissueSiteGrid S : dataWriters){
             S.writeData(writer,gridNo);
             gridNo++;
        }

    }
    public void Close(){
        if(summary!=null) summary.Close();

    }

}
