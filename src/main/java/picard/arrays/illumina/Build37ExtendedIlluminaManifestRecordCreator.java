package picard.arrays.illumina;

import htsjdk.samtools.reference.ReferenceSequenceFile;
import picard.PicardException;

import java.io.File;
import java.util.Map;

public class Build37ExtendedIlluminaManifestRecordCreator {

    private final Map<String, ReferenceSequenceFile> referenceFilesMap;
    private final Map<String, File> chainFilesMap;

    Build37ExtendedIlluminaManifestRecordCreator(final Map<String, ReferenceSequenceFile> referenceFilesMap,
                                                 final Map<String, File> chainFilesMap) {
        this.referenceFilesMap = referenceFilesMap;
        this.chainFilesMap = chainFilesMap;
    }

    public Build37ExtendedIlluminaManifestRecord validateLocusEntryAndCreateExtendedRecord(final IlluminaBPMLocusEntry locusEntry,
                                                                                       final IlluminaManifestRecord record) {
        if (!locusEntry.ilmnId.equals(record.getIlmnId())) {
            throw new PicardException("Mismatch");
        }
        return new Build37ExtendedIlluminaManifestRecord(record, Build37ExtendedIlluminaManifestRecord.Flag.PASS, "1", 0, "A", "C", "T", "rs22");
    }

}
