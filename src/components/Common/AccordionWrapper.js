import React, { useState } from 'react';
import PropTypes from 'prop-types';
import {
  AccordionTitle,
  ContentAndButtonWrapper,
  ContentWrapper,
  CollapseButton,
  Title,
  HelpButton,
  HelpText,
} from './Accordion.styles';

const AccordionWrapper = ({
  children,
  maxHeight,
  title,
  defaultCollapsed,
  defaultShowHelp,
  helpButtonText,
  helpText,
  horizontal = false,
}) => {
  const [state, setState] = useState({
    collapsed: defaultCollapsed,
    showHelp: defaultShowHelp,
  });

  const toggleCollapsed = (e) => {
    e.preventDefault();
    setState({ collapsed: !state.collapsed });
  };

  const toggleHelp = (e) => {
    e.preventDefault();
    setState({ showHelp: !state.showHelp });
  };

  const style = horizontal ? { height: '100%' } : { width: '100%' };

  return (
    <div style={style}>
      <Title>
        <CollapseButton onClick={toggleCollapsed}>
          {state.collapsed ? '+' : '-'}
        </CollapseButton>
        <AccordionTitle>{title}</AccordionTitle>
        {helpText !== null && (
          <HelpButton onClick={toggleHelp}>{helpButtonText}</HelpButton>
        )}
      </Title>
      <ContentAndButtonWrapper>
        <ContentWrapper maxHeight={maxHeight} collapsed={state.collapsed}>
          {helpText !== null && (
            <HelpText show={state.showHelp}>{helpText}</HelpText>
          )}
          {children}
        </ContentWrapper>
      </ContentAndButtonWrapper>
    </div>
  );
};

AccordionWrapper.propTypes = {
  children: PropTypes.oneOfType([
    PropTypes.arrayOf(PropTypes.node),
    PropTypes.node,
  ]),
  maxHeight: PropTypes.string.isRequired,
  title: PropTypes.node.isRequired,
  defaultCollapsed: PropTypes.bool,
  defaultShowHelp: PropTypes.bool,
  helpButtonText: PropTypes.string,
  helpText: PropTypes.oneOfType([
    PropTypes.arrayOf(PropTypes.node),
    PropTypes.node,
  ]),
};
AccordionWrapper.defaultProps = {
  defaultCollapsed: false,
  defaultShowHelp: false,
  helpButtonText: 'Show Help',
  helpText: null,
};

export default AccordionWrapper;
