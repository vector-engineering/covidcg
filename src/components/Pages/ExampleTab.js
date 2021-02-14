import React from 'react';
import { observer } from 'mobx-react';
import styled from 'styled-components';
import { useStores } from '../../stores/connect';
import _ from 'underscore';

import { getGene } from '../../utils/gene_protein';
import { getLocationByNameAndLevel } from '../../utils/location';
import { queryPrimers } from '../../utils/primer';
import { ISOToInt } from '../../utils/date';

import {
  GROUP_SNV,
  DNA_OR_AA,
  COORDINATE_MODES,
  NORM_MODES,
  COUNT_MODES,
  DATE_BINS,
  ASYNC_STATES,
  TABS,
} from '../../constants/defs.json';
import { config } from '../../config';

import ExternalLink from '../Common/ExternalLink';
import SkeletonElement from '../Common/SkeletonElement';
// import LoadingSpinner from '../Common/LoadingSpinner';
import TempImage from '../../assets/images/cg_short_v13@4x_square.png';

import GlobalLineagesImage from '../../assets/analysis_screens/global_lineages.png';
// import XinfadiImage from '../../assets/analysis_screens/xinfadi.png';
import IcelandImage from '../../assets/analysis_screens/iceland.png';
import D614GWestCoastImage from '../../assets/analysis_screens/d614g_west_coast.png';
import D614GUSStatesImage from '../../assets/analysis_screens/d614g_us_states.png';
import USCDCPrimerImage from '../../assets/analysis_screens/us_cdc_primer.png';
import D614GEuropeNAImage from '../../assets/analysis_screens/d614g_europe_na.png';
import N203204Image from '../../assets/analysis_screens/n_203_204_coocurrence.png';

import ImageExample1_1 from '../../assets/example_screens/example1_1.png';
import ImageExample1_2 from '../../assets/example_screens/example1_2.png';
import ImageExample1_3 from '../../assets/example_screens/example1_3.png';
import ImageExample1_4 from '../../assets/example_screens/example1_4.png';

import ImageExample2_1 from '../../assets/example_screens/example2_1.png';
import ImageExample2_2 from '../../assets/example_screens/example2_2.png';
import ImageExample2_3 from '../../assets/example_screens/example2_3.png';
//import ImageExample2_4 from '../../assets/example_screens/example2_4.png';

import ImageExample3_1 from '../../assets/example_screens/example3_1.png';
import ImageExample3_2 from '../../assets/example_screens/example3_2.png';
import ImageExample3_3 from '../../assets/example_screens/example3_3.png';

const ExampleTabContainer = styled.div`
  background-color: #f8f8f8;
`;

const ExampleTabContent = styled.div`
  display: flex;
  flex-direction: column;
  align-items: stretch;

  flex-grow: 1;
  max-width: 1000px;
  padding: 20px;
  margin: 0 auto;

  background-color: #fff;
`;

const ExampleHeader = styled.div`
  padding-left: 10px;
  max-width: 800px;

  p {
    font-weight: normal;
    line-height: normal;
  }
`;

const ExampleTitle = styled.h2``;

const ExampleList = styled.div`
  display: flex;
  flex-direction: row;
  flex-wrap: wrap;
`;
const ExampleItem = styled.a`
  display: flex;
  flex-direction: column;
  align-items: stretch;

  width: 300px;
  height: 250px;
  margin: 10px;
  border: 1px solid #ccc;
  box-shadow: 0px 3px 4px #aaa;

  text-decoration: none;

  transition: 0.1s all ease-in-out;
  &:hover {
    border: 1px solid #00e;
  }
`;
const ExampleItemImage = styled.div`
  display: flex;
  flex-direction: row;
  align-items: center;
  justify-content: center;
  overflow: hidden;

  img {
    width: auto;
    height: 150px;
  }
`;
const ExampleItemFooter = styled.div`
  height: 100px;
  padding: 10px;

  border-top: 1px solid #aaa;

  .example-item-title {
    display: block;
    font-size: 1.1em;
    font-weight: 700;
    color: #444;
  }

  .example-item-description {
    display: block;
    font-size: 0.85em;
    font-weight: normal;
    color: #888;
    margin: 3px 0px;
    line-height: normal;
  }
`;

const TOC = styled.div`
  font-weight: normal;
`;

const ExampleTutorial = styled.div`
  padding: 10px;
  padding-top: 0px;
  font-weight: normal;
  line-height: normal;
  max-width: 800px;
`;

const TutorialImage = styled.img`
  width: 750px;
  border: 1px solid #aaa;
`;

const scrollToRef = (id, e) => {
  if (e !== undefined) {
    e.preventDefault();
  }
  const el = window.document.getElementById(id);
  window.scrollTo(0, el.offsetTop);
};

const ExampleTab = observer(() => {
  const {
    configStore,
    plotSettingsStore,
    UIStore,
    locationDataStore,
  } = useStores();

  const selectTree = locationDataStore.selectTree;
  const exampleItems = [
    {
      title: 'Global Lineages',
      description: 'View the growth of the B lineage family over all locations',
      image: GlobalLineagesImage,
      settings: {
        plotSettings: {
          groupStackNormMode: NORM_MODES.NORM_PERCENTAGES,
          groupStackCountMode: COUNT_MODES.COUNT_NEW,
          groupStackDateBin: DATE_BINS.DATE_BIN_WEEK,
        },
        UI: {
          activeTab: TABS.TAB_GROUP,
        },
        config: {
          groupKey: config.group_cols.lineage.name,
          dnaOrAa: DNA_OR_AA.DNA,
          selectedGene: getGene('S'),
          coordinateMode: COORDINATE_MODES.COORD_GENE,
          selectedLocationNodes: [selectTree], // select root
        },
      },
    },
    // {
    //   title: 'Lineages in China - Beijing Xinfadi Market',
    //   description:
    //     "New lineages in uncovered in early June in Beijing's Xinfadi Market may have been circulating in China in March",
    //   image: XinfadiImage,
    //   settings: {
    //     plotSettings: {
    //       groupStackNormMode: NORM_MODES.NORM_COUNTS,
    //       groupStackCountMode: COUNT_MODES.COUNT_NEW,
    //       groupStackDateBin: DATE_BINS.DATE_BIN_DAY,
    //     },
    //     UI: {
    //       activeTab: TABS.TAB_GROUP,
    //     },
    //     config: {
    //       groupKey: config.group_cols.lineage.name,
    //       dnaOrAa: DNA_OR_AA.AA,
    //       selectedGene: getGene('S'),
    //       coordinateMode: COORDINATE_MODES.COORD_GENE,
    //       selectedLocationNodes: [getLocationByNameAndLevel(selectTree, 'China', 'country')[0]],
    //     },
    //   },
    // },
    {
      title: 'Lineages in Iceland',
      description:
        'Lineages sequenced in Iceland during the early stages of the pandemic',
      image: IcelandImage,
      settings: {
        plotSettings: {
          groupStackNormMode: NORM_MODES.NORM_COUNTS,
          groupStackCountMode: COUNT_MODES.COUNT_NEW,
          groupStackDateBin: DATE_BINS.DATE_BIN_DAY,
        },
        UI: {
          activeTab: TABS.TAB_GROUP,
        },
        config: {
          groupKey: config.group_cols.lineage.name,
          dnaOrAa: DNA_OR_AA.AA,
          selectedGene: getGene('S'),
          coordinateMode: COORDINATE_MODES.COORD_GENE,
          selectedLocationNodes: [
            getLocationByNameAndLevel(selectTree, 'Iceland', 'country')[0],
          ],
        },
      },
    },
    {
      title: 'Rise of Spike D614G mutation in West Coast USA',
      description:
        'The Spike D614G mutation has accumulated in frequency in the West Coast states of the USA',
      image: D614GWestCoastImage,
      settings: {
        plotSettings: {
          groupStackNormMode: NORM_MODES.NORM_PERCENTAGES,
          groupStackCountMode: COUNT_MODES.COUNT_NEW,
          groupStackDateBin: DATE_BINS.DATE_BIN_DAY,
        },
        UI: {
          activeTab: TABS.TAB_GROUP,
        },
        config: {
          groupKey: GROUP_SNV,
          dnaOrAa: DNA_OR_AA.AA,
          selectedGene: getGene('S'),
          coordinateMode: COORDINATE_MODES.COORD_GENE,
          selectedLocationNodes: [
            getLocationByNameAndLevel(selectTree, 'Washington', 'division')[0],
            getLocationByNameAndLevel(selectTree, 'Oregon', 'division')[0],
            getLocationByNameAndLevel(selectTree, 'California', 'division')[0],
          ],
          dateRange: [ISOToInt('2020-03-01'), ISOToInt('2020-06-01')],
          selectedGroups: [{ group: 'S|614|D|G' }],
        },
      },
    },
    {
      title: 'Prevalence of Spike D614G in various US States',
      description:
        'The proportion of sequences with the Spike D614G mutation varies between US States',
      image: D614GUSStatesImage,
      settings: {
        plotSettings: {
          locationDateNormMode: NORM_MODES.NORM_PERCENTAGES,
          locationDateCountMode: COUNT_MODES.COUNT_CUMULATIVE,
          locationDateDateBin: DATE_BINS.DATE_BIN_DAY,
        },
        UI: {
          activeTab: TABS.TAB_LOCATION,
        },
        config: {
          groupKey: GROUP_SNV,
          dnaOrAa: DNA_OR_AA.AA,
          selectedGene: getGene('S'),
          coordinateMode: COORDINATE_MODES.COORD_GENE,
          selectedGroups: [{ group: 'S|614|D|G' }],
          selectedLocationNodes: [
            getLocationByNameAndLevel(selectTree, 'Washington', 'division')[0],
            getLocationByNameAndLevel(selectTree, 'California', 'division')[0],
            getLocationByNameAndLevel(selectTree, 'Florida', 'division')[0],
            getLocationByNameAndLevel(
              selectTree,
              'Massachusetts',
              'division'
            )[0],
            getLocationByNameAndLevel(selectTree, 'Michigan', 'division')[0],
            getLocationByNameAndLevel(selectTree, 'New York', 'division')[0],
            getLocationByNameAndLevel(selectTree, 'Texas', 'division')[0],
            getLocationByNameAndLevel(selectTree, 'Wisconsin', 'division')[0],
          ],
        },
      },
    },
    {
      title: 'Diagnostics: NT SNVs in US CDC qPCR primer/probe sequences',
      description:
        'Prevalence of any mutations present within the US CDC primer and probe sequences (N1 + N2), for sequences in the US',
      image: USCDCPrimerImage,
      settings: {
        plotSettings: {
          groupStackNormMode: NORM_MODES.NORM_PERCENTAGES,
          groupStackCountMode: COUNT_MODES.COUNT_CUMULATIVE,
          groupStackDateBin: DATE_BINS.DATE_BIN_DAY,
        },
        UI: {
          activeTab: TABS.TAB_GROUP,
        },
        config: {
          groupKey: GROUP_SNV,
          dnaOrAa: DNA_OR_AA.DNA,
          coordinateMode: COORDINATE_MODES.COORD_PRIMER,
          selectedLocationNodes: [
            getLocationByNameAndLevel(selectTree, 'USA', 'country')[0],
          ],
          selectedPrimers: []
            .concat(
              queryPrimers({ Institution: 'US CDC', Name: '2019-nCoV-N1-F' })
            )
            .concat(
              queryPrimers({ Institution: 'US CDC', Name: '2019-nCoV-N1-R' })
            )
            .concat(
              queryPrimers({ Institution: 'US CDC', Name: '2019-nCoV-N1-P' })
            )
            .concat(
              queryPrimers({ Institution: 'US CDC', Name: '2019-nCoV-N2-F' })
            )
            .concat(
              queryPrimers({ Institution: 'US CDC', Name: '2019-nCoV-N2-R' })
            )
            .concat(
              queryPrimers({ Institution: 'US CDC', Name: '2019-nCoV-N2-P' })
            ),
        },
      },
    },
    {
      title: 'Co-occurrence of R203K and G204R in N gene',
      description:
        'Two SNVs in the N gene, R203K and G204R, co-occur with each other. These SNVs are associated with the B.1.1 lineage.',
      image: N203204Image,
      settings: {
        plotSettings: {
          groupStackNormMode: NORM_MODES.NORM_PERCENTAGES,
          groupStackCountMode: COUNT_MODES.COUNT_NEW,
          groupStackDateBin: DATE_BINS.DATE_BIN_DAY,
        },
        UI: {
          activeTab: TABS.TAB_GROUP,
        },
        config: {
          groupKey: GROUP_SNV,
          dnaOrAa: DNA_OR_AA.AA,
          selectedGene: getGene('N'),
          coordinateMode: COORDINATE_MODES.COORD_GENE,
          selectedLocationNodes: [
            getLocationByNameAndLevel(selectTree, 'USA', 'country')[0],
          ],
          selectedGroups: [{ group: 'N|203|R|K' }, { group: 'N|204|G|R' }],
        },
      },
    },
    {
      title: 'Emergence of Spike D614G in Europe/North America',
      description: '',
      image: D614GEuropeNAImage,
      settings: {
        plotSettings: {
          groupStackNormMode: NORM_MODES.NORM_PERCENTAGES,
          groupStackCountMode: COUNT_MODES.COUNT_NEW,
          groupStackDateBin: DATE_BINS.DATE_BIN_WEEK,
        },
        UI: {
          activeTab: TABS.TAB_GROUP,
        },
        config: {
          groupKey: GROUP_SNV,
          dnaOrAa: DNA_OR_AA.AA,
          selectedGene: getGene('S'),
          coordinateMode: COORDINATE_MODES.COORD_GENE,
          selectedLocationNodes: [
            getLocationByNameAndLevel(selectTree, 'North America', 'region')[0],
            getLocationByNameAndLevel(selectTree, 'Europe', 'region')[0],
          ],
          selectedGroups: [{ group: 'S|614|D|G' }],
        },
      },
    },
  ];

  const onExampleClick = (title, e) => {
    e.preventDefault();

    const exampleItem = _.findWhere(exampleItems, { title });
    // console.log(exampleItem);

    // Apply settings for each store
    Object.keys(exampleItem.settings).forEach((store) => {
      // Call the example action for each store
      if (store === 'config') {
        configStore.resetValues(exampleItem.settings[store]);
      } else if (store === 'plotSettings') {
        plotSettingsStore.resetValues(exampleItem.settings[store]);
      } else if (store === 'UI') {
        UIStore.resetValues(exampleItem.settings[store]);
      }
    });
  };

  const exampleElements = [];
  exampleItems.forEach((exampleItem) => {
    exampleElements.push(
      <ExampleItem
        key={`example-item-${exampleItem.title}`}
        href="#"
        onClick={onExampleClick.bind(this, exampleItem.title)}
      >
        <ExampleItemImage>
          <img
            src={
              exampleItem.image === undefined ? TempImage : exampleItem.image
            }
          />
        </ExampleItemImage>
        <ExampleItemFooter>
          <span className="example-item-title">{exampleItem.title}</span>
          <p className="example-item-description">{exampleItem.description}</p>
        </ExampleItemFooter>
      </ExampleItem>
    );
  });

  const renderExamples = () => {
    // Hide the examples while the app is still initializing
    if (UIStore.caseDataState === ASYNC_STATES.STARTED) {
      const skeletonList = [];
      for (let i = 0; i < exampleElements.length; i++) {
        skeletonList.push(
          <SkeletonElement
            key={`example-loading-${i}`}
            delay={2}
            width="300px"
            height={250}
            style={{
              margin: '10px',
            }}
          />
        );
      }
      return skeletonList;
    }

    return exampleElements;
  };

  return (
    <ExampleTabContainer>
      <a id="getting-started-top" />
      <ExampleTabContent>
        <TOC>
          <h3>Table of Contents</h3>
          <ul>
            <li>
              <a
                href="#introduction"
                onClick={scrollToRef.bind(this, 'introduction')}
              >
                <b>Introduction</b>
              </a>
            </li>
            <li>
              <a
                href="#tutorial-1"
                onClick={scrollToRef.bind(this, 'tutorial-1')}
              >
                <b>Tutorial</b>: Tracking the new S477N mutation in Australia
              </a>
            </li>
            <li>
              <a
                href="#tutorial-2"
                onClick={scrollToRef.bind(this, 'tutorial-2')}
              >
                <b>Tutorial</b>: Checking diagnostic primers or probes for SNVs
              </a>
            </li>
            <li>
              <a
                href="#tutorial-3"
                onClick={scrollToRef.bind(this, 'tutorial-3')}
              >
                <b>Tutorial</b>: Tracking the Spike D614G mutant distribution in
                your location of interest over time
              </a>
            </li>
            <li>
              <a
                href="#example-analyses"
                onClick={scrollToRef.bind(this, 'example-analyses')}
              >
                <b>Example Analyses</b>
              </a>
            </li>
          </ul>
        </TOC>

        <ExampleHeader style={{ marginBottom: '10px' }}>
          <a id="introduction" />
          <ExampleTitle>Introduction</ExampleTitle>
          <a href="#" onClick={scrollToRef.bind(this, 'getting-started-top')}>
            Back
          </a>
        </ExampleHeader>

        <ExampleTutorial>
          <b style={{ fontSize: '1.4em', marginTop: '10px' }}>
            Preprint:{' '}
            <ExternalLink href="https://www.biorxiv.org/content/10.1101/2020.09.23.310565v2">
              https://www.biorxiv.org/content/10.1101/2020.09.23.310565v2
            </ExternalLink>
          </b>
          <p>
            <b>
              The COVID-19 CoV Genetics browser was designed to empower diverse
              projects on hCoV-19 (SARS-CoV-2) transmission, evolution,
              emergence, immune interactions, diagnostics, therapeutics,
              vaccines, and tracking of interventions.
            </b>
          </p>
          <p>
            Tracking the evolution of the emerging coronavirus is essential for
            scientists and public health professionals, as well as developers of
            vaccines, diagnostics, and therapeutics. This work is enabled by
            data generously shared from contributors across the world via{' '}
            <ExternalLink href="https://www.gisaid.org/">GISAID</ExternalLink>,
            on which this research is based.
          </p>
          <p>
            COVID-19 CG helps users to quickly find answers to questions
            including but not limited to:
          </p>
          <ol>
            <li>
              Which clades and lineages are present in a given city or region
              within a user-specified period of time?
            </li>
            <li>
              Which variants should I test my therapeutic, antibody, or
              diagnostic on before implementation in a specific region?
            </li>
            <li>
              What are the community outcomes after particular policies,
              vaccines, or therapeutics are applied in that population?
            </li>
            <li>
              Are there data from transient mutations that can elucidate common
              mechanisms of resistance to acquired immunity? Can this be
              leveraged for vaccine, antibody, or small molecule drug design?
            </li>
          </ol>
          <p>
            Users can view the comprehensive nucleotide and amino acid residue
            variation in their selection to inform their research hypothesis
            generation or anti-COVID-19 product development. For example,
            COVID-19 CG enables users to evaluate commonly used or custom
            primers/probes or targets/epitopes based on their location and dates
            of interest.
            {/* As sequencing efforts start to include details about hCoV-19 (SARS-CoV-2)
          isolates, users can also sort virus data according to patient
          characteristics such as age, gender, clinical status, isolate type, as
          well as passaging, sequencing, and assembly method. */}
          </p>
          <p>
            To accelerate COVID-19 research and public health efforts, COVID-19
            CG will be continually upgraded with new and improved features so
            that users can quickly and reliably pinpoint critical mutations as
            the virus evolves throughout the pandemic.
          </p>
          <p>
            Towards this goal, we strongly advocate that countries continue to
            generate data sampled from patients, animals, and environments and
            share their data in a timely manner via{' '}
            <ExternalLink href="https://www.gisaid.org/">GISAID</ExternalLink>,
            so that scientists across the world are maximally informed about
            developments in the spread of the emerging coronavirus responsible
            for COVID-19.
          </p>
          <p>
            Reach out to us{' '}
            <ExternalLink href="https://twitter.com/covidcg">
              @covidcg
            </ExternalLink>{' '}
            on twitter, or email us at{' '}
            <ExternalLink href="mailto:covidcg@broadinstitute.org">
              covidcg@broadinstitute.org
            </ExternalLink>
          </p>
        </ExampleTutorial>

        <ExampleHeader>
          <a id="tutorial-1" />
          <ExampleTitle>
            Tutorial: Tracking the new S477N mutation in Australia
          </ExampleTitle>
          <a href="#" onClick={scrollToRef.bind(this, 'getting-started-top')}>
            Back
          </a>
        </ExampleHeader>

        <ExampleTutorial>
          <p>
            In this example, we will show you how to track mutations of interest
            in your locations of interest. There is a new dominant SARS-CoV-2
            variant in Australia as you can see in the COVID-19 CG{' '}
            <b>Compare Lineages</b> tab.
          </p>
          <TutorialImage src={ImageExample1_1} />
          <p>
            To select the <b>B.1.1.25</b> lineage, click on it in the{' '}
            <b>Legend</b>, or click on the bars representing the lineage in the{' '}
            <b>Lineage Plot</b>.
          </p>
          <TutorialImage src={ImageExample1_2} />
          <p>
            Now that we know there could be a unique mutation in the Spike,
            S477N, that sets apart this lineage, B.1.1.25, let’s check it out in
            the <b>Compare AA SNVs</b> tab! You can get to this tab by simply
            switching from Lineage to SNV view in the navigational menu (top
            left sidebar).
          </p>
          <TutorialImage src={ImageExample1_3} />
          <p>
            Let’s find out how many sequences bearing the S477N mutation have
            been detected in different parts of Australia in the{' '}
            <b>Compare Locations</b> tab. One tip is to deselect one location
            with the least number of SARS-CoV-2 sequences - this tells the
            COVID-19 CG site that you want to compare locations within
            Australia. From the data, it looks like the vast majority of
            S477N-harboring sequences were obtained from Victoria, Australia
            (and the earliest sequence was obtained in January!).
          </p>
          <TutorialImage src={ImageExample1_4} />
          <p>
            Alternatively, pick other countries to compare to Australia. For
            S477N, we recommend the United Kingdom and Florida. Although, please
            note that COVID-19 CG only reflects data contributed to GISAID.
            Variants of interest could be present in other countries, and at
            earlier dates, but not yet known to the public because the
            sequencing centers in those countries have not collected or
            deposited their data in GISAID.
          </p>
        </ExampleTutorial>

        <ExampleHeader>
          <a id="tutorial-2" />
          <ExampleTitle>
            Tutorial: Checking diagnostic primers or probes for SNVs
          </ExampleTitle>
          <a href="#" onClick={scrollToRef.bind(this, 'getting-started-top')}>
            Back
          </a>
        </ExampleHeader>

        <ExampleTutorial>
          <p>
            In this example, we will show you how to check your diagnostic
            primer or probe for SNVs present in your targeted location of
            diagnostic implementation. Using this approach, we found at least 11
            diagnostic primer pairs commonly used around the world that could be
            impacted by SNVs found near their 3’ end.
          </p>
          <TutorialImage src={ImageExample2_1} />
          <TutorialImage src={ImageExample2_2} />
          <p>
            To figure out where these SNVs are located globally, you can click
            on the SNV of interest (e.g., G10097A) and switch to the{' '}
            <b>Compare Locations</b> tab - simply click on the{' '}
            <b>Compare Locations</b> tab! For G10097A, we found, using COVID-19
            CG, that it has been detected in every continent but is particularly
            common in some countries.
          </p>
          <TutorialImage src={ImageExample2_3} />
          <p>
            This image was downloaded in the <b>Compare Locations</b> tab.
          </p>
          <p>
            We advocate that labs and clinics use COVID-19 CG (
            <ExternalLink href="https://covidcg.org">
              https://covidcg.org
            </ExternalLink>
            ) to check their most commonly used primers and probes against the
            SARS-CoV-2 sequences that are prevalent in their geographic regions.
            More than 700 commonly used primers or probes are now available for
            COVID-19 CG users to select under the <b>Primers/Probes</b> menu.
          </p>
        </ExampleTutorial>

        <ExampleHeader>
          <a id="tutorial-3" />
          <ExampleTitle>
            Tutorial: Tracking the Spike D614G mutant distribution in your
            location of interest over time
          </ExampleTitle>
          <a href="#" onClick={scrollToRef.bind(this, 'getting-started-top')}>
            Back
          </a>
        </ExampleHeader>

        <ExampleTutorial>
          <p>
            In this example, we will show you how to track the distribution of
            the D614G mutation across SARS-CoV-2 isolates in your locations of
            interest over time.
          </p>
          <p>
            Warning: we strongly caution against inferring (i) chains or
            directionality of transmission and (ii) changes in the
            transmissibility of any SARS-CoV-2 SNV based on population dynamics
            alone. Inconsistent sampling, sampling biases, differences in
            founder host population traits (even median patient age),
            superspreading events, regionally and temporally differential travel
            restrictions, and numerous other factors instead of virus biological
            differences can influence the global distribution of SNVs.
          </p>
          <p>
            Let’s take a look at the city-island-country of Singapore.
            SARS-CoV-2 sequences carrying the D614G mutation are shown in pink!
            In other words, any variant that does not carry the G614 AA SNV is
            shown in grey.
          </p>
          <TutorialImage src={ImageExample3_1} />
          <p>
            Let’s change things up and take a look at South Korea. The
            population dynamics of the D614G SNV can look dramatically different
            between locations. In Singapore, the D614 SNV appeared more common
            over time as the bulk of the G614 variants waned by April. In South
            Korea, the opposite trend is apparent, with the bulk of the D614
            variants subsiding by April, and a wave of G614 variants appearing
            in May through July.
          </p>
          <p>
            We encourage you to explore other countries or regions on your own.
          </p>
          <TutorialImage src={ImageExample3_2} />
          <p>
            How about comparing the percent of SARS-CoV-2 sequences that carry
            the G614 mutation by locations of interest? Here, we picked three
            Indian cities in the <b>Compare Locations</b> tab. The{' '}
            <b>Location-Date Plot</b> can be adjusted to show Cumulative or New
            sequences, Percentages or Counts, grouped by Day, Week, or Month.
          </p>
          <TutorialImage src={ImageExample3_3} />
        </ExampleTutorial>

        <ExampleHeader>
          <a id="example-analyses" />
          <ExampleTitle>Example Analyses</ExampleTitle>
          <a href="#" onClick={scrollToRef.bind(this, 'getting-started-top')}>
            Back
          </a>
          <p>
            Use these example analyses to get started and explore the features
            of this application. If you would like to add an analysis to this
            list, please{' '}
            <ExternalLink href="https://github.com/vector-engineering/covidcg">
              submit a pull request on our GitHub
            </ExternalLink>
            , or contact us at{' '}
            <ExternalLink href="mailto:covidcg@broadinstitute.org">
              covidcg@broadinstitute.org
            </ExternalLink>
          </p>
        </ExampleHeader>
        <ExampleList>{renderExamples()}</ExampleList>
      </ExampleTabContent>
    </ExampleTabContainer>
  );
});

export default ExampleTab;
