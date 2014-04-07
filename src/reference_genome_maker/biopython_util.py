"""
Methods that encapsulate genome editing techniques applied to BioPython
objects.

NOTE: These functions are copied from a bigger set of utils currently in
private development in the Church Lab. Eventually we'll migrate to depending
on a common library.
"""

from Bio.SeqFeature import FeatureLocation
from Bio.SeqFeature import SeqFeature


def add_feature_to_seq_record(seq_record, new_seq_feature, ignore_ids=False):
    """Adds the feature to the SeqRecord list of features and updates
    the custom mapping that we use.

    Args:
        seq_record: The SeqRecord to modify.
        new_seq_feature: The SeqFeature to add.
        ignore_ids: If True, does not check for unique ids. Only works when
            not using an underlying map.

    Raises:
        AssertionError if a feature with that id already exists.
    """
    if (hasattr(seq_record, 'feature_id_to_index_map') and ignore_ids):
        raise AssertionError(
                "Cannot ignore ids when using feature_id_to_index_map in "
                "SeqRecord.")

    # HACK: Assuming our feature ids are safe, then this feature does not
    # need to be added again.
    if (hasattr(seq_record, 'feature_id_to_index_map') and
            new_seq_feature.id in seq_record.feature_id_to_index_map):
        return

    not_unique_err = "Feature id already exists %s" % (new_seq_feature.id,)
    if hasattr(seq_record, 'feature_id_to_index_map'):
        assert new_seq_feature.id not in seq_record.feature_id_to_index_map,\
                not_unique_err
        seq_record.features.append(new_seq_feature)
        last_index = len(seq_record.features) - 1
        seq_record.feature_id_to_index_map[new_seq_feature.id] = last_index
        assert len(seq_record.feature_id_to_index_map.keys()) == len(seq_record.features)
    else:
        if ignore_ids:
            seq_record.features.append(new_seq_feature)
            return

        # Otherwise check uniqueness.
        unique = True
        for feature in seq_record.features:
            if feature.id == new_seq_feature.id:
                unique = False
                break
        if unique:
            seq_record.features.append(new_seq_feature)
        else:
            raise AssertionError(not_unique_err)


def delete_interval(seq_record, interval, validation_seq=None):
    """Delete all of the bases in the interval and update features
    appropriately. Mutates seq_record.

    When updating features:
        The feature's end position is shifted upstream for each deleted
        base that overlapped the feature.

    Args:
        seq_record: SeqRecord object that will be mutated to reflect changes.
        interval: Pythonic interval to delete (inclusive lower-bound, exclusive
            upper-bound). Assummes a 0-indexed genome.
        validation_seq: If provided, we validate that this is the sequence
            in the interval. For debugging purposes.

    Returns:
        The mutated SeqRecord.
    """
    if validation_seq:
        assert validation_seq == str(
                seq_record.seq[interval[0]:interval[1]]), (
                "\tValidation sequence doesn't match interval.\n"
                "\tInterval: %s\n"
                "\tActual seq (truncated): %s\n"
                "\tActual seq start +/- 5: %s\n"
                "\tValidation seq: %s" % (
                        str(interval),
                        str(seq_record.seq[interval[0]:
                                min(interval[0] + 6, interval[1])]),
                        str(seq_record.seq[interval[0] - 5:interval[0]]) + '*' +
                                str(seq_record.seq[interval[0]:
                                        interval[0] + 5]),
                        validation_seq
                )
        )

    interval_size = interval[1] - interval[0]

    # Change the underlying sequence.
    new_seq = (
            seq_record.seq[:interval[0]] +
            seq_record.seq[interval[0] + interval_size:]
    )
    seq_record.seq = new_seq

    # Update the features.
    updated_features = []
    for feature in seq_record.features:
        if feature.location.start >= interval[1]:
            # Entire feature is downstream of deletion. Move the whole thing.
            feature = feature._shift(-1 * interval_size)
        elif does_interval_overlap_feature(interval, feature):
            # Feature overlaps interval.
            feature = _shift_feature_overlapped_by_event(feature, interval, -1)
        if feature:
            updated_features.append(feature)
    seq_record.features = updated_features

    # Return the mutated seq_record.
    return seq_record


def does_interval_overlap_feature(interval, feature):
    """Checks whether the given interval overlaps the feature's location.

    Args:
        interval: A two-tuple of integers (start, end).
        feature: A SeqFeature.

    Returns:
        A boolean indicating whether the features overlap.
    """
    interval_start = interval[0]
    interval_end = interval[1]

    if feature.location.start == interval_start:
        # >>>>>>>>>
        # (....)
        return interval_end - interval_start > 0

    elif feature.location.start < interval_start:
        # >>>>>>>>>
        #    (..........)
        return feature.location.end > interval_start

    else:
        #      >>>>>>>>>
        # (........)
        return feature.location.start < interval_end


def _shift_feature_overlapped_by_event(feature, event_interval, event_polarity):
    """Shifts the feature and its sub_features depending on which parts
    of the interval overlap it.

    NOTE: Clients should use the returned feature.

    Args:
        feature: SeqFeature that may be mutated, but clients should
            use the return value.
        event_interval: The interval of the event in the frame of the
            SeqRecord to which the feature belongs.
        event_polarity: 1 for insertion, -1 for deletion.

    Returns:
        An updated copy of the feature, or None.
    """
    assert event_polarity == -1 or event_polarity == 1

    event_interval_size = event_interval[1] - event_interval[0]

    # Determine how much to shift the feature
    if (event_interval[0] <= feature.location.start and
            event_interval[1] >= feature.location.end):
        # Event covers entire feature. Remove it.
        #   >>>>>>>>>
        # (...........)
        return None
    elif (event_interval[0] > feature.location.start and
            event_interval[1] < feature.location.end):
        # Event is inside feature. We'll shift the end
        # >>>>>>>>>>
        #   (....)
        shift_amount = event_interval_size
    elif (event_interval[0] > feature.location.start and
            event_interval[1] >= feature.location.end):
        # >>>>>>>>>
        #   (........)
        shift_amount = feature.location.end - event_interval[0]
    else:
        #     >>>>>>>>>
        # (........)
        shift_amount = event_interval[1] - feature.location.start
    _shift_end_of_feature(feature, shift_amount * event_polarity)

    # Update the sub_features appropriately.
    updated_sub_features = []
    for sub_feature in feature.sub_features:
        if sub_feature.location.start >= event_interval[1]:
            # Entire feature is downstream of deletion. Move the whole thing.
            sub_feature = sub_feature._shift(event_polarity * event_interval_size)
        elif does_interval_overlap_feature(event_interval, sub_feature):
            # Feature overlaps interval.
            sub_feature = _shift_feature_overlapped_by_event(
                    sub_feature, event_interval, event_polarity)
        if sub_feature:
            updated_sub_features.append(sub_feature)
    feature.sub_features = updated_sub_features

    return feature


def _shift_end_of_feature(feature, shift_amount):
    """Mutates the feature by shifting the end of its location.
    """
    feature.location = FeatureLocation(
            start=feature.location.start,
            end=feature.location.end._shift(shift_amount),
            strand=feature.strand)


def insert_sequence_and_update_features(seq_record, insert_seq,
        insert_start_position, extend_feature_ends=False,
        no_feature_shift_if_inside=False,
        insert_feature_type=None,
        insert_feature_id=None,
        insert_feature_strand=1):
    """Inserts a sequence at the given position and updates feature locations.

    This method differs from insert_sequence() in that this method doesn't
    "know" what it's being used for.

    Features are updated according to the following rules:
        * If a feature lies upstream of an insertion, nothing to do.
        * If the insertion is inside of a feature, increase the size of the
            feature the appropriate amount.
            * If the feature has sub_features, then recursively apply
                these rules.

    Args:
        seq_record: The mutable SeqRecord object.
        insert_seq: The sequence to insert, as read along the forward strand.
        insert_start_position: The start position of the insertion. Pythonic.
        extend_feature_ends: If True, any features whose end position is at
            insert_start_position will have their end extended to capture
            the entire insertion.
        no_feature_shift_if_inside: Don't shift the feature even if the change
            is inside the feature.

    Returns:
        The mutated seq_record.
    """
    # Update the sequence.
    insert_size = len(insert_seq)
    new_seq = (
            seq_record.seq[:insert_start_position] +
            insert_seq +
            seq_record.seq[insert_start_position:]
    )
    assert len(new_seq) == len(seq_record.seq) + insert_size
    seq_record.seq = new_seq

    # Update the features.
    updated_features = []
    for feature in seq_record.features:
        if feature.location.start > insert_start_position:
            # Feature is downstream of deletion. Move the whole thing.
            feature = feature._shift(insert_size)
        elif no_feature_shift_if_inside:
            pass
        elif ((feature.location.end > insert_start_position and
                    feature.location.start <= insert_start_position) or
                (extend_feature_ends and
                        feature.location.end == insert_start_position)):
            # Insertion is inside the feature.
            _shift_feature_with_event_at_position(
                    feature, insert_start_position, insert_size,
                    extend_feature_ends=extend_feature_ends)
        updated_features.append(feature)
    seq_record.features = updated_features

    # If insert_feature_type is defined, add a new feature to represent the inserted
    # sequence.
    if insert_feature_type and insert_feature_id:
        new_seq_feature = SeqFeature(
                location=FeatureLocation(insert_start_position,
                        insert_start_position + insert_size),
                type=insert_feature_type,
                strand=insert_feature_strand,
                id=insert_feature_id
        )
        new_seq_feature.qualifiers[
                FeatureQualifierExtensionKeys.GENOME_REFACTOR_ID] = (
                        insert_feature_id)
        add_feature_to_seq_record(seq_record, new_seq_feature)


def _shift_feature_with_event_at_position(
        feature, event_position, shift_amount, extend_feature_ends=False):
    """Mutates the feature by shifting its end the given amount.

    The Bio.SeqFeature._shift() method doesn't update the sub-features
    correctly which is probably why it is underscore-indicated-private.

    So the effect of this method is:
        1) For each sub-feature:
            * Shift the whole thing, if after the event.
            * Shift just the end, if overlapping the event.
            * No shift, if the sub_feature lies entirely before the event.
        2) For the whole feature:
            * Shift the end of the whole feature.
    """
    # Save original length for final validation.
    # NOTE: When a SeqFeature has sub_features, SeqFeature.__len__() returns
    # the sum of the lengths of the sub-features rather than the size of
    # the feature.location interval, so we are explicit here with the assertion
    # check.
    orig_feature_len = len(feature.location)

    ### 1) Update the sub_features appropriately.
    updated_sub_features = []
    for sub_feature in feature.sub_features:
        if sub_feature.location.start > event_position:
            sub_feature = sub_feature._shift(shift_amount)
        elif ((sub_feature.location.end > event_position and
                    sub_feature.location.start <= event_position) or
                (extend_feature_ends and
                        sub_feature.location.end == event_position)):
            _shift_end_of_feature(sub_feature, shift_amount)
        updated_sub_features.append(sub_feature)
    feature.sub_features = updated_sub_features

    ### 2) Shift the end of the feature.
    _shift_end_of_feature(feature, shift_amount)

    # Assert the operation succeeded.
    expected_updated_feature_len = orig_feature_len + shift_amount
    assert expected_updated_feature_len == len(feature.location), (
            "Feature: %s, Expected: %d, actual: %d" % (
                    feature,
                    expected_updated_feature_len,
                    len(feature)))
