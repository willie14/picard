package picard.arrays.illumina;

import htsjdk.samtools.reference.ReferenceSequenceFile;
import picard.PicardException;

import java.io.File;
import java.util.Map;

public class ExtendedIlluminaManifestRecordCreator {

    private final Map<String, ReferenceSequenceFile> referenceFilesMap;
    private final Map<String, File> chainFilesMap;

    ExtendedIlluminaManifestRecordCreator(final Map<String, ReferenceSequenceFile> referenceFilesMap,
                                          final Map<String, File> chainFilesMap) {
        this.referenceFilesMap = referenceFilesMap;
        this.chainFilesMap = chainFilesMap;
    }

    public ExtendedIlluminaManifestRecord validateLocusEntryAndCreateExtendedRecord(final IlluminaBPMLocusEntry locusEntry,
                                                                                    final IlluminaManifestRecord illuminaManifestRecord,
                                                                                    final boolean dupe) {
        ExtendedIlluminaManifestRecord record = new ExtendedIlluminaManifestRecord(locusEntry, illuminaManifestRecord, dupe);
        return record;
    }

}
