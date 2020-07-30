import styled from 'styled-components';

const EmptyPlot = styled.div`
  display: flex;
  flex-direction: row;
  align-items: center;
  justify-content: center;

  height: ${({ height }) => height}px;

  font-size: 1.5em;
  font-weight: normal;
  color: #888;
  line-height: normal;

  padding: 10px 30px;

  p {
    text-align: center;
    max-width: 800px;
  }
`;
EmptyPlot.defaultProps = {
  height: 250,
};

export default EmptyPlot;
