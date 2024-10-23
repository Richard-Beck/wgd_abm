package Examples._3OffLatticeExample;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;

import HAL.Gui.OpenGL2DWindow;
import HAL.Tools.FileIO;
import java.io.*;

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
        outputPath = outputHome.getAbsolutePath() + "/run_" + folderList.length;
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
        summary.Write("time,normal,tumor,wgd,vessel,immune\n");
    }

    IO(){
        disable_output = true;
    };

    public void RecordOut(ExampleOffLattice E, int time){
        if(disable_output) return;
        int normal = 0, tumor=0, wgd= 0, vessel=0, immune=0;
        for (CellOL cell : E) {
            if(cell.type==CellOL.CellType.tumor) tumor++;
            if(cell.type==CellOL.CellType.immune) immune++;
            if(cell.type==CellOL.CellType.vessel) vessel++;
            if(cell.type==CellOL.CellType.normal) normal++;
            if(cell.type==CellOL.CellType.wgd) wgd++;
        }
        summary.Write(time+","+normal+","+tumor+","+wgd+","+
                vessel+","+immune+"\n");
    }

    public void Draw(int time, OpenGL2DWindow vis){
        if(disable_output) return;
        if(time%imageWriteFrequency!=0) return;
        if(vis==null) return;
        vis.ToJPG(imageFolder.getAbsolutePath()+"/"+time+".jpg");
    }

    public void Close(){
        if(summary!=null) summary.Close();

    }

}
