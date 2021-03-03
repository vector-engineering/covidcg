import styled from 'styled-components';

export const FooterContainer = styled.div`
  margin-top: auto;
  display: flex;
  flex-direction: column;
  align-items: stretch;
  background-color: #f8f8f8;
  flex-shrink: 0;

  padding: 5px 8px;
  border-top: 1px solid #ccc;

  font-size: 0.7rem;
  font-weight: normal;
  line-height: normal;
`;

export const DAA = styled.div`
  margin-bottom: 5px;
  padding-bottom: 5px;
  border-bottom: 1px solid #ccc;
`;

export const Version = styled.div`
  margin-bottom: 2px;
  span.version-num {
    font-weight: bold;
  }
  a {
    margin-left: 3px;
  }
`;

export const SequenceMeta = styled.div`
  span.date {
    font-weight: bold;
  }
`;
