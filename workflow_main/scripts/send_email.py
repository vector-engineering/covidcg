#!/usr/bin/env python3
# coding: utf-8

"""using SendGrid's Python Library
https://github.com/sendgrid/sendgrid-python

Author: Albert Chen (Deverman Lab, Broad Institute)
"""

import argparse
import base64
import os
import sys

from functools import reduce
from pathlib import Path
from sendgrid import SendGridAPIClient
from sendgrid.helpers.mail import (
    From,
    To,
    Mail,
    Attachment,
    FileContent,
    FileName,
    FileType,
    Disposition,
    ContentId,
)


def send_email(recipients, subject, html_content, attachments=None):
    """Send an email via. SendGrid

    Parameters
    ----------
    recipients: str or list of str
    subject: str
    html_content: str
    attachments: (Optional) list of str
        List of file paths

    Returns
    -------
    None
    """

    # Expand recipients, if they were passed in as one comma-delimited string
    recipients = reduce(lambda x, y: x + y.split(","), recipients, [])

    message = Mail(
        from_email="vector-engineering@broadinstitute.org",
        to_emails=recipients,
        subject=subject,
        html_content=html_content,
    )

    message_attachments = []
    for file_path in attachments:

        with open(file_path, "rb") as f:
            data = f.read()
            f.close()

        encoded = base64.b64encode(data).decode()
        attachment = Attachment()
        attachment.file_content = FileContent(encoded)
        attachment.file_type = FileType("text/plain")
        attachment.file_name = FileName(Path(file_path).name)
        attachment.disposition = Disposition("attachment")
        # attachment.content_id = ContentId('Example Content ID')

        message_attachments.append(attachment)

    if message_attachments:
        message.attachment = message_attachments

    try:
        sg = SendGridAPIClient(os.environ.get("SENDGRID_API_KEY"))
        response = sg.send(message)
        print(response.status_code)
        print(response.body)
        print(response.headers)
    except Exception as e:
        print(e.message)
        sys.exit(1)
        pass


def main():
    """Command-line entry point"""

    # Load command line arguments
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--recipients",
        type=str,
        nargs="+",
        required=True,
        default=[],
        help="Email recipients",
    )
    parser.add_argument("--subject", type=str, required=True, help="Email Subject")
    parser.add_argument(
        "--content", type=str, required=True, help="Email content (HTML)"
    )
    parser.add_argument(
        "--attachments", type=str, nargs="+", default=[], help="Path to attachments"
    )

    args = parser.parse_args()

    send_email(args.recipients, args.subject, args.content, args.attachments)


if __name__ == "__main__":
    main()
