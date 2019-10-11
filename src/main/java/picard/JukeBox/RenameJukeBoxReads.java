package picard.JukeBox;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.reference.ReferenceSequence;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.PicardException;
import picard.analysis.SinglePassSamProgram;
import picard.cmdline.programgroups.BaseCallingProgramGroup;
import picard.util.CsvInputParser;

import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;

@CommandLineProgramProperties(

        summary = "nope",
        oneLineSummary = "nope",
        programGroup = BaseCallingProgramGroup.class
)

public class RenameJukeBoxReads extends SinglePassSamProgram {
    @Argument
    public String RUN_BARCODE;
    @Argument
    public File TILES_DESCRIPTION;
    @Argument
    public File BEADS_XY;

    private CsvInputParser csvParser;
    private DataInputStream cooredinateReader;
    private JukeBoxTile currentTileRow;
    private SAMFileWriter outputWriter;

    @Override
    protected void setup(final SAMFileHeader header, final File samFile) {
        csvParser = new CsvInputParser(false, TILES_DESCRIPTION);
        try {
            cooredinateReader = new DataInputStream(new FileInputStream(BEADS_XY));
        } catch (FileNotFoundException e) {
            throw new PicardException("Trouble opening file: " + BEADS_XY, e);
        }
        final String[] headerRow = csvParser.next();
        if (!headerRow[0].equals("TileID") ||
                !headerRow[1].equals("TileOffset") ||
                !headerRow[2].equals("BeadsPerTile") ||
                !headerRow[3].equals("Ring") ||
                !headerRow[4].equals("Radius") ||
                !headerRow[5].equals("Theta")) {
            throw new PicardException("Unexpected format in header row of file: " + TILES_DESCRIPTION);
        }
        UpdateTileRow();
        final SAMFileHeader newHeader = header.clone();
        newHeader.setSortOrder(SAMFileHeader.SortOrder.unsorted);
        outputWriter = new SAMFileWriterFactory().makeWriter(newHeader, true, OUTPUT, referenceSequence.getReferenceFile());
    }

    @Override
    protected void acceptRead(final SAMRecord rec, final ReferenceSequence ref) {
        final String indexString = rec.getReadName();
        final long index = Long.parseLong(indexString);
        while (index >= currentTileRow.TileOffset + currentTileRow.BeadsPerTile) {
            UpdateTileRow();
        }
        final JukeBoxCoordinates readCoordinates;
        try {
            readCoordinates = getNextCoordinates();
        } catch (IOException e) {
            throw new PicardException("Trouble reading from binary file: " + BEADS_XY, e);
        }
        final String newReadname = String.format("%s:%s:%d:%d:%d",
                RUN_BARCODE, indexString, currentTileRow.TileID,
                readCoordinates.xCoordinate, readCoordinates.yCoordinate);

        rec.setReadName(newReadname);
        outputWriter.addAlignment(rec);
    }

    private void UpdateTileRow() {
        currentTileRow = new JukeBoxTile(csvParser.next());
    }

    private final byte[] buffer = new byte[4];

    private JukeBoxCoordinates getNextCoordinates() throws IOException {
        final float x = getLittleEndianFloat();
        final float y = getLittleEndianFloat();

        return new JukeBoxCoordinates(x, y);
    }

    private float getLittleEndianFloat() throws IOException {
        if(cooredinateReader.read(buffer, 0, 4) !=4){
            throw new PicardException("Unable to read 4 bytes from binary file: " + BEADS_XY );
        }
        return ByteBuffer.wrap(buffer).order(ByteOrder.LITTLE_ENDIAN).getFloat();
    }

    @Override
    protected void finish() {
        csvParser.close();
        outputWriter.close();

        try {
            cooredinateReader.close();
        } catch (IOException e) {
            throw new PicardException("Trouble closing file: " + BEADS_XY, e);
        }
    }

    static private class JukeBoxTile {
        final long TileID;
        final long TileOffset;
        final long BeadsPerTile;
        final long Ring;
        final long Radius;
        final float Theta;

        JukeBoxTile(String[] row) {
            this.TileID = Long.parseLong(row[0]);
            this.TileOffset = Long.parseLong(row[1]);
            this.BeadsPerTile = Long.parseLong(row[2]);
            this.Ring = Long.parseLong(row[3]);
            this.Radius = Long.parseLong(row[4]);
            this.Theta = Float.parseFloat(row[5]);
        }
    }

    static private class JukeBoxCoordinates {
        final int xCoordinate;
        final int yCoordinate;
        private static final int MULTIPLIER = 10;

        JukeBoxCoordinates(float x, float y) {
            this.xCoordinate = Math.round(x * MULTIPLIER);
            this.yCoordinate = Math.round(y * MULTIPLIER);
        }
    }
}
