import React from 'react';
import PropTypes from 'prop-types';

import {
  WalkthroughListContainer,
  WalkthroughListTitle,
  WalkthroughItemList,
  WalkthroughItemContainer,
  WalkthroughItemImage,
  WalkthroughItemText,
  WalkthroughItemTitle,
  WalkthroughItemDescription,
} from './WalkthroughList.styles';

const WalkthroughItem = ({ title, description, image, link }) => {
  return (
    <WalkthroughItemContainer
      title={title}
      href={link}
      target="_blank"
      rel="noopener noreferrer"
    >
      <WalkthroughItemImage image={image} />
      <WalkthroughItemText>
        <WalkthroughItemTitle>{title}</WalkthroughItemTitle>
        <WalkthroughItemDescription>{description}</WalkthroughItemDescription>
      </WalkthroughItemText>
    </WalkthroughItemContainer>
  );
};
WalkthroughItem.propTypes = {
  title: PropTypes.string.isRequired,
  description: PropTypes.string.isRequired,
  image: PropTypes.string.isRequired,
  link: PropTypes.string.isRequired,
};

const walkthroughs = [
  {
    title: 'COVID CG Intro',
    description: 'Introduction to the site',
    image:
      'https://storage.googleapis.com/ve-public/walkthrough/20220119_CG_intro_banner.jpg',
    link: 'https://storage.googleapis.com/ve-public/walkthrough/20220119_CG_intro.pdf',
  },
  {
    title: 'Case study: Paxlovid + Mpro (nsp5)',
    description:
      'Do prominent variants have mutations that might impact the Pfizer protease inhibitor?',
    image:
      'https://storage.googleapis.com/ve-public/walkthrough/20220119_CG_case_study_mpro_banner.jpg',
    link: 'https://storage.googleapis.com/ve-public/walkthrough/20220119_CG_case_study_mpro.pdf',
  },
  {
    title: 'Case study: Ab binding on Omicron spike',
    description:
      'Could mutations in the Omicron spike RBD impact therapeutic antibody binding?',
    image:
      'https://storage.googleapis.com/ve-public/walkthrough/20220119_CG_case_study_omicron_spike_banner.jpg',
    link: 'https://storage.googleapis.com/ve-public/walkthrough/20220119_CG_case_study_omicron_spike.pdf',
  },
  {
    title: 'Case study: Mutations in CDC primers',
    description:
      'Could any mutations be affecting the CDC N diagnostic primers?',
    image:
      'https://storage.googleapis.com/ve-public/walkthrough/20220119_CG_case_study_CDCN1_banner.jpg',
    link: 'https://storage.googleapis.com/ve-public/walkthrough/20220119_CG_case_study_CDCN1.pdf',
  },
];

const WalkthroughList = () => {
  const walkthroughItems = [];

  walkthroughs.forEach((walkthrough, i) => {
    walkthroughItems.push(
      <WalkthroughItem key={`walkthrough-item-${i}`} {...walkthrough} />
    );
  });

  return (
    <WalkthroughListContainer>
      <WalkthroughListTitle>Walkthroughs</WalkthroughListTitle>
      <WalkthroughItemList>{walkthroughItems}</WalkthroughItemList>
    </WalkthroughListContainer>
  );
};

export default WalkthroughList;
