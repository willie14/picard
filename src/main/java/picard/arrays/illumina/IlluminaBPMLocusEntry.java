package picard.arrays.illumina;

import htsjdk.tribble.annotation.Strand;

/**
 * A simple class to represent a locus entry in an Illumina Bead Pool Manifest (BPM) file
 */
public class IlluminaBPMLocusEntry {
    int version;        // The LocusEntry version.

    // IlmnID (probe identifier) of locus
    String ilmnId;      // Verified.

    // Name (variant identifier) of locus
    String name;        // Verified

    // Index of this entry.
    int index;

    // Illumina Strand value
    IlluminaManifestRecord.IlluminaStrand ilmnStrand;       // Verified

    // SNP value for locus (e.g., [A/C])
    String snp;     // Verified

    // Chromosome for the locus (e.g., XY)
    String chrom;   // V

    String ploidy;  // V

    String species; // V

    // Mapping location of locus
    int mapInfo;    // V

    // Customer Strand
    String customerStrand;      // ?? Not verified - NOT in csv?

    // AddressA ID of locus
    int addressA;       // V

    // Only populated in CSV files or BPM files with version 4 data block
    String alleleAProbeSeq;     // V

    // AddressB ID of locus (0 if none)
    int addressB;       // V

    // Only populated in CSV files or BPM files with version 4 data block (empty if none)
    String alleleBProbeSeq;     // V

    String genomeBuild;         // V
    String source;              // V
    String sourceVersion;       // V
    IlluminaManifestRecord.IlluminaStrand sourceStrand;     // V

    // Only populated in CSV files or BPM files with version 4 data block
    String sourceSeq;           // V

    // Only populated in CSV files or BPM files with version 4 data block
    String topGenomicSeq;       // V

    int expClusters;            // V
    boolean intensityOnly;          // V

   // Identifies type of assay (0 - Infinium II , 1 - Infinium I (A/T), 2 - Infinium I (G/C)
    int assayType;              // Not Verified - NOT in csv

    float fracA;              // Not Verified - NOT in csv
    float fracC;              // Not Verified - NOT in csv
    float fracT;              // Not Verified - NOT in csv
    float fracG;              // Not Verified - NOT in csv

    // Refstrand annotation
    Strand refStrand;           // V

    // Not part of the locusEntry record in the BPM, added here for convenience
    int normalizationId;        // Not Verified - added here.

    public IlluminaBPMLocusEntry() {
        version = -1;

        ilmnId = "";
        name = "";
        index = -1;
        ilmnStrand = IlluminaManifestRecord.IlluminaStrand.NONE;
        snp = "";
        chrom = "";
        ploidy = "";
        species = "";
        mapInfo = -1;
        customerStrand = "";

        addressA = -1;
        addressB = -1;

        genomeBuild = "";
        source = "";
        sourceVersion = "";
        sourceStrand = IlluminaManifestRecord.IlluminaStrand.NONE;

        expClusters = -1;
        intensityOnly = false;
        assayType = -1;

        fracA = 0.0f;
        fracC = 0.0f;
        fracT = 0.0f;
        fracG = 0.0f;

        refStrand = Strand.NONE;

        normalizationId = -1;
    }

}
