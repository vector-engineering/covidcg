/*
 * Send a request to the Google Cloud Function to decrypt our hashed Accession IDs
 */

export function decryptAccessionIds(accessionIds) {
  return fetch(
    'https://us-central1-test-project-1-2-3.cloudfunctions.net/decrypt-accession-id',
    {
      method: 'POST',
      mode: 'cors',
      cache: 'no-cache',
      credentials: 'same-origin',
      headers: { 'Content-Type': 'application/json' },
      redirect: 'follow',
      referrerPolicy: 'no-referrer-when-downgrade',
      body: JSON.stringify({
        accession_ids: accessionIds,
      }),
    }
  )
    .then((response) => response.json())
    .catch((error) => console.error(error));
}
