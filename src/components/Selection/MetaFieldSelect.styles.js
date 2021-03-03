import styled from 'styled-components';

const formWidth = '160px';

export const MetaFieldSelectContainer = styled.div`
  display: flex;
  flex-direction: column;
  align-items: stretch;

  padding-left: 15px;
  margin-bottom: 10px;

  span.title {
    font-weight: 500;
    font-size: 1rem;
    margin-bottom: 5px;
  }
`;

export const SelectList = styled.div`
  padding-left: 10px;
`;

export const SelectContainer = styled.form`
  display: flex;
  flex-direction: row;
  align-items: center;
  justify-content: flex-start;

  margin-top: 3px;
  padding-right: 15px;

  label {
    margin-right: 5px;
    font-weight: normal;
    width: 11em;
  }

  .metadata-multi-select {
    // flex-grow: 1;
    width: ${formWidth};

    .dropdown-heading {
      height: 20px;
      padding: 0 8px;

      .dropdown-heading-value {
        font-weight: normal;
        .gray {
          color: #666;
        }
      }
    }

    .dropdown-content {
      width: 250px;
      right: 0;
      padding-top: 0px;

      z-index: 2;
      .panel-content {
        border-radius: 0px;
        .select-item {
          padding: 3px 8px;
          font-weight: normal;
        }
      }
    }
  }
`;

// const PatientAgeContainer = styled.div`
//   display: flex;
//   flex-direction: row;
//   align-items: center;

//   width: ${formWidth};
//   // padding-left: 5px;

//   span {
//     font-weight: normal;
//   }
//   span.from-label {
//   }
//   span.to-label {
//     margin-left: 5px;
//   }
//   input {
//     margin-left: 5px;
//   }
// `;
