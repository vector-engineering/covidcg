import { getAckTextsFromAccessionIds } from './acknowledgements';

// https://riptutorial.com/javascript/example/24711/client-side-csv-download-using-blob
// https://stackoverflow.com/questions/14964035/how-to-export-javascript-array-info-to-csv-on-client-side
// function downloadCsv(csvString, filename) {
//   let blob = new Blob([csvString]);
//   let link = window.document.createElement('a');
//   let url = URL.createObjectURL(blob);
//   link.setAttribute('href', url);
//   link.setAttribute('download', filename);
//   link.style.visibility = 'hidden';
//   parent.window.document.body.appendChild(link);
//   link.click();
//   parent.window.document.body.removeChild(link);
// }

function downloadAcknowledgements(selectedAccessionIds) {
  // Get the list of selected Accession IDs, and map to
  // acknowledgement texts
  let ackTexts = getAckTextsFromAccessionIds(selectedAccessionIds);
  // console.log(ackTexts);

  // Write to a CSV string
  // Accession ID and sample date first
  // Then acknowledgement texts
  let csvString =
    'Accession ID,Collection Date,Originating lab,Submitting lab,Authors\n';

  for (let i = 0; i < selectedAccessionIds.length; i++) {
    // Write Accession ID
    csvString += selectedAccessionIds[i]['gisaid_id'] + ',';
    // Write Sample Date
    // Get the date in ISO format, and chop off the time/timezone info at the end
    // So that we get YYYY-MM-DD (the same as the original input format)
    csvString +=
      new Date(selectedAccessionIds[i]['sample_date'])
        .toISOString()
        .substring(0, 10) + ',';

    // Write Acknowledgement texts
    // Since these can contain commas, wrap each in double quotes
    csvString += '"' + ackTexts[i]['originating_lab'] + '",';
    csvString += '"' + ackTexts[i]['submitting_lab'] + '",';
    csvString += '"' + ackTexts[i]['authors'] + '"\n';
  }

  //let a = new WorkerGlobalScope();
  //console.log(self);

  //downloadCsv(csvString, 'acknowledgements.csv');

  let blob = new Blob([csvString]);
  let url = URL.createObjectURL(blob);

  return {
    blobURL: url,
  };
}

self.addEventListener(
  'message',
  function (e) {
    const data = JSON.parse(e.data);
    //console.log('in downloadworker event listener', data);

    let result = 'asdf';
    if (data.type === 'downloadAcknowledgements') {
      // This is a terminal endpoint, we don't need to post a message back
      result = downloadAcknowledgements(data.selectedAccessionIds);
    }
    // console.log(result);
    self.postMessage(JSON.stringify(result));
  },
  false
);
