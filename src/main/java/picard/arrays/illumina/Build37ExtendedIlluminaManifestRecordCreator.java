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
        validateEntryField(locusEntry.ilmnId, record.getIlmnId(), "ilmnId");
        validateEntryField(locusEntry.name, record.getName(), "name");
        validateEntryField(locusEntry.ilmnStrand, record.getIlmnStrand(), "ilmnStrand");
        validateEntryField(locusEntry.snp, record.getSnp(), "snp");
        validateEntryField(locusEntry.chrom, record.getChr(), "chrom");
        validateEntryField(locusEntry.ploidy, record.getPloidy(), "ploidy");
        validateEntryField(locusEntry.species, record.getSpecies(), "species");
        validateEntryField(locusEntry.mapInfo, record.getPosition(), "mapInfo");
        validateEntryField(locusEntry.addressA, Integer.parseInt(record.getAddressAId()), "addressAId");
        if (locusEntry.version == 4) {
            validateEntryField(locusEntry.alleleAProbeSeq, record.getAlleleAProbeSeq(), "alleleAProbeSeq");
        }
        if ((locusEntry.addressB != -1) && (record.getAddressBId() != null)) {
            validateEntryField(locusEntry.addressB, Integer.parseInt(record.getAddressBId()), "addressBId");
        }
        if (locusEntry.version == 4) {
            validateEntryField(locusEntry.alleleBProbeSeq, record.getAlleleBProbeSeq(), "alleleBProbeSeq");
        }
        validateEntryField(locusEntry.genomeBuild, record.getGenomeBuild(), "genomeBuild");
        validateEntryField(locusEntry.source, record.getSource(), "source");
        validateEntryField(locusEntry.sourceVersion, record.getSourceVersion(), "sourceVersion");
        validateEntryField(locusEntry.sourceStrand, record.getSourceStrand(), "sourceStrand");
        if (locusEntry.version == 4) {
            validateEntryField(locusEntry.sourceSeq, record.getSourceSeq(), "sourceSeq");
            validateEntryField(locusEntry.topGenomicSeq, record.getTopGenomicSeq(), "topGenomicSeq");
        }
        if (record.getExpClusters() != null) {
            validateEntryField(locusEntry.expClusters, Integer.parseInt(record.getExpClusters()), "expClusters");
        }
        validateEntryField(locusEntry.intensityOnly, record.getIntensityOnly(), "intensityOnly");
        if (locusEntry.version == 8) {
            validateEntryField(locusEntry.refStrand, record.getRefStrand(), "refStrand");
        }

        return new Build37ExtendedIlluminaManifestRecord(record, Build37ExtendedIlluminaManifestRecord.Flag.PASS, "1", 0, "A", "C", "T", "rs22");
    }

    private void validateEntryField(final Object locusEntryField, final Object recordField, final String fieldName) {
        if (!locusEntryField.equals(recordField)) {
            throw new PicardException("Field '" + fieldName + "' disagrees between BPM file (found '" + locusEntryField + "') and CSV (found: '" + recordField + "')");
        }
    }

}
