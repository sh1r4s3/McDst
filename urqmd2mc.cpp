/*
 * vim:et:sw=4:
 *
 * Copyright (c) eternity Nikita (sh1r4s3) Ermakov <sh1r4s3@mail.si-head.nl>
 *
 * SPDX-License-Identifier: Beerware
 *
 *
 * urqmd2mc reads UrQMD events from the ftn13 or ftn14 ascii files and
 * converts them to the UniGen format and saves on a root file.
 *
 * ftn14 contains snapshots at given times (event steps). The event steps
 * are converted to separate events.
 *
 * ftn13 contains the final snapshot and the freeze-out coordinates.
 * The freeze-out coordinates are used. The final snapshot coordinates
 * are discarded.
 */

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <map>
#include <algorithm>
#include <optional>
#include <stdexcept>

#include "TObject.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TBranch.h"
#include "TMath.h"
#include "TClonesArray.h"

#include "McRun.h"
#include "McEvent.h"
#include "McParticle.h"
#include "McPIDConverter.h"
#include "McArrays.h"

TFile *fi;
TTree *tr;
TClonesArray *arrays[McArrays::NAllMcArrays];

void bomb(const char *myst) {
    std::cerr << "Error: " << myst << ", bombing" << std::endl;
    exit(-1);
}


std::string newName(char* origName) {
    std::string fname(origName);
    std::string key1 = ".f13";
    std::string key2 = ".f14";
    std::size_t found1 = fname.rfind(key1);
    std::size_t found2 = fname.rfind(key2);
    if ( found1 != std::string::npos ) {
        fname.replace(found1, key1.length(), ".mcDst.root");
    }
    else if ( found2 != std::string::npos ) {
        fname.replace(found2, key2.length(), ".mcDst.root");
    }
    else {
        bomb("Wrong input data format (not f13 or f14)");
    }

    return fname;
}

int trapco(int ityp, int ichg) {
    // translate UrQMD pid code to pdg code

    /* UrQMD PIDs are in fact composite - a particle is fully defined by the
       type specifier (ityp), the charge (ichg) and in case of baryons, the
       third component isospin (iso3; ignored by us at present). For
       simplicity, our conversion tables collapse these into a single number
       as follows:
       - propagate the sign of ityp (particle-antiparticle distinction for
       baryons, strangeness-antistrangeness distinction for mesons) to that
       of the stored value;
       - shift the ichg range from -2..2 (UrQMD does not support nuclear
       fragments other than protons and neutrons so all particles it
       produces fall in this range) to 0..4 to make sure it doesn't
       interfere with the above;
       - multiply shifted charge by +/-1000 and add it to type. The latter
       is guaranteed to be smaller than 1000 (baryon types are one- or
       two-digit, meson types three-digit) so that way no ambiguities
       occur. */
    int id;
    if (ityp >= 0) id = 1000 * (ichg + 2) + ityp;
    else
        id = -1000 * (ichg + 2) + ityp;

    /* This will return 0 for unknown input values. */
    return McPIDConverter::instance()->pdgCode(id, McPIDConverter::eUrQMD);
}

bool optExists(int argc, char ** argv, const std::string & opt) {
    return std::find(argv, argv + argc, opt) != argv + argc;
}

std::optional<std::string> getOptArg(int argc, char ** argv, const std::string & opt) {
    auto itArg = std::find(argv, argv + argc, opt);
    if (itArg == argv + argc || ++itArg == argv + argc) {
        return std::nullopt;
    }

    return std::make_optional(*itArg);
}

void usage() {
    std::cerr << "Usage: urqmd2mc -i <input file> -o <output file> [-v]" << std::endl
              << "       optional arguments:" << std::endl
              << "       -v  verbose mode" << std::endl;
}

int main(int argc, char ** argv) {
    std::ifstream in;
    char c;
    std::string dust;

    // Switcher that excludes elastic collisions
    bool excludeElastic = false;

    // McRun initialization
    McRun *run = nullptr;

    bool isElastic;
    std::string version, comment;
    int filetype, eos, aproj, zproj, atarg, ztarg, nr;
    double beta, b, bmin, bmax, sigma, elab, plab, sqrts, time, dtime;

    if (argc < 5) { // for input, output and number of events
        usage();
        return 1;
    }

    auto inputFile = getOptArg(argc, argv, "-i");
    auto outputFile = getOptArg(argc, argv, "-o");

    if (!inputFile || !outputFile) { // we require input and output file
        usage();
        return 1;
    }

    bool verbose = optExists(argc, argv, "-v");

    if (verbose) {
        std::cout << "input file: " << *inputFile << std::endl
                  << "output file: " << *outputFile << std::endl;
    }

    unsigned nout {0};

    // Check that filename contains .f13 or .f14
    TString oFileName(*outputFile);

    // Try to open file
    in.open(*inputFile);
    if (in.fail()) {
        bomb("cannot open input file");
    }

    fi = TFile::Open(oFileName.Data(), "RECREATE", "UrQMD");
    fi->SetCompressionLevel(9);
    int bufsize = 65536 * 4;
    int split = 99;

    // Create and set out TTree
    McEvent::Class()->IgnoreTObjectStreamer();
    McParticle::Class()->IgnoreTObjectStreamer();
    McRun::Class()->IgnoreTObjectStreamer();

    tr = new TTree("McDst", "UrQMD tree", split);
    tr->SetAutoSave(1000000);
    for (unsigned int i = 0; i < McArrays::NAllMcArrays; i++) {
        // Create arrayss
        arrays[i] = new TClonesArray( McArrays::mcArrayTypes[i], McArrays::mcArraySizes[i] );
        // Set branch
        tr->Branch(McArrays::mcArrayNames[i], &arrays[i], bufsize / 4, split);
    }

    int events_processed {0};
    unsigned ntime_slices {0};

    // Start event loop
    while (!in.eof()) {
        std::string line;

        // Read event information
        in >> dust >> dust >> version >> dust >> dust;
        in >> dust >> filetype >> dust >> dust >> dust >> aproj >> zproj;
        in >> dust >> dust >> dust >> atarg >> ztarg;
        in >> dust >> dust >> dust >> beta >> dust >> dust;
        in >> dust >> b >> bmin >> bmax >> dust >> sigma;
        in >> dust >> eos >> dust >> elab >> dust >> sqrts >> dust >> plab;
        in >> dust >> nr >> dust >> dust >> dust;
        in >> dust >> dust >> time >> dust >> dtime;
        in.ignore(777,'\n'); // ignore the rest of the line

        if (verbose) {
            std::cout << "event#: " << events_processed
                      << "version: " << version
                      << "sqrts: " << sqrts
                      << "dtime: " << dtime << std::endl;
        }

        comment.clear();
        // read 4 lines of options and 6 lines of params
        for (int i=0; i<10; i++) {
            getline(in,line);
            comment.append(line);
            comment.append("\n");
        }
        in.ignore(777,'\n');

        // Clear all arrays
        for (unsigned int i = 0; i < McArrays::NAllMcArrays; i++) {
            // Create arrayss
            arrays[i]->Clear();
        }

        // Increment number of processed events
        events_processed++;

        int step_nr=0;
        char pee;
        // Loop over time slices
        while (!in.eof()) {
            int mult;
            double step_time;
            pee=in.peek();
            if (pee=='U') break;
            if (pee==EOF) break;
            in >> mult >> step_time;

            // Check if collision is elastic
            isElastic = false;
            if ( (aproj+atarg) == mult ) {
                isElastic = true;
            }

            if ( verbose ) {
                std::cout << "Number of particles in event: " << mult << std::endl;
            }

            in.ignore(777,'\n'); // ignore the rest of the line
            getline(in,line);

            // Loop over generated particles
            for (int i=0;  i<mult; i++) {
                if ( verbose ) {
                    std::cout << "Working on particle i: " << i << std::endl;
                }
                double t, x, y, z, e, px, py, pz;
                int ityp, iso3, ichg, status, parent, parent_decay, mate;
                int decay, child[2];

                // Read particl information
                in >> t >> x >> y >> z;
                in >> e >> px >> py >> pz >> dust;
                in >> ityp >> iso3 >> ichg >> mate >> dust >> dust;

                // Print particle information
                if ( verbose ) {
                    std::cout << Form( " t: %6.3f \tx: %6.3f \ty: %6.3f \tz: %6.3f \tpx: %6.3f \tpy: %6.3f \tpz: %6.3f\n",
                                      t, x, y, z, px, py, pz );
                }

                // Read freeze-out information
                if (filetype==13) {
                    in >> t >> x >> y >> z;
                    in >> e >> px >> py >> pz;

                    // Print freeze-out information
                    if ( verbose ) {
                        std::cout << Form( " t: %6.3f \tx: %6.3f \ty: %6.3f \tz: %6.3f \tpx: %6.3f \tpy: %6.3f \tpz: %6.3f\n",
                                          t, x, y, z, px, py, pz );
                    }
                } // if (filetype==13)

                if (in.fail()) bomb("while reading tracks");

                // Do not fill McParticles for elastic collisions (if skip)
                if (isElastic && excludeElastic) continue;

                status = parent_decay = decay = child[0] = child[1] = 0;
                // Add new particle to the event
                McParticle * k = new McParticle(i, trapco(ityp, ichg), status, parent,
                                                       parent_decay, mate-1, decay, child,
                                                       px, py, pz, e, x, y, z, t );
                arrays[McArrays::Particle]->operator[](i) = k;

                // Print particle information stored in McParticle
                if (verbose) {
                    int iPart = arrays[McArrays::Particle]->GetEntries();
                    std::cout << "PART=" << iPart << std::endl;
                    McParticle *particle = (McParticle*)arrays[McArrays::Particle]->At(iPart-1);
                    if ( !particle ) {
                        std::cout << "Particle does not exist!" << std::endl;
                        break;
                    }
                    std::cout << Form( " t: %6.3f \tx: %6.3f \ty: %6.3f \tz: %6.3f \tpx: %6.3f \tpy: %6.3f \tpz: %6.3f\n",
                                      particle->t(), particle->x(),
                                      particle->y(), particle->z(),
                                      particle->px(), particle->py(),
                                      particle->pz() );
                }
            }

            do in.get(c); while (c!='\n');

            // Create an instance of the McEvent and add it to the DST
            // only when the event is not elastic and those are not required
            // to be skipped
            if ( !(isElastic && excludeElastic) ) {
                // Add new McEvent
                McEvent * ev = new McEvent();
                ev->setEventNr(nr);
                ev->setB(b);
                ev->setPhi(0);
                ev->setNes((int) (time/dtime));
                ev->setComment(line.data());
                ev->setStepNr(step_nr++);
                ev->setStepT(step_time);

                std::cout << "event#: " << events_processed << "    tslice: " << ntime_slices << std::endl;
                int iEvent = arrays[McArrays::Event]->GetEntries();
                arrays[McArrays::Particle]->operator[](iEvent) = ev;
                ntime_slices++;
                // Fill DST with event and track information
                nout += tr->Fill();
            }
        }

        if (pee==EOF) break;
    }
    in.close();
    std::cout << events_processed << " events processed\n";

    tr->Write();

    // create the run object
    std::string generator = "UrQMD";
    generator.append(version);
    double m = 0.938271998;
    double ecm = sqrts/2; // energy per nucleon in cm
    double pcm = sqrt(ecm*ecm-m*m); // momentum per nucleon in cm
    double gamma = 1.0/sqrt(1-beta*beta);
    double pproj = gamma*(+pcm-beta*ecm);
    double ptarg = gamma*(-pcm-beta*ecm);
    run = new McRun( generator.data(), comment.data(),
                    aproj, zproj, pproj,
                    atarg, ztarg, ptarg,
                    bmin, bmax, -1, 0, 0, sigma, events_processed);
    run->Write();
    fi->Write();
    fi->Close();
    std::cout << "Total bytes were written: " << nout << std::endl;
    return 0;
}
